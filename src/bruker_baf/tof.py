import sqlite3
import struct
from contextlib import closing
from dataclasses import dataclass
from math import isfinite, sqrt
from pathlib import Path


@dataclass(frozen=True)
class TofCalibrationRecord:
    start_time: float
    time_step: float
    time_offset: float
    quadratic_coefficient: float
    mass_time_correction: float
    higher_order_correction: float
    mass_offset: float
    raw_blob: bytes

    @property
    def mode(self) -> int:
        return 10

    @classmethod
    def from_bytes(cls, data: bytes) -> "TofCalibrationRecord":
        if len(data) < 4:
            raise ValueError("Incomplete TOF calibration blob")
        mode = struct.unpack_from("<i", data)[0]
        if mode != 10:
            raise NotImplementedError(
                f"Unsupported TOF calibration mode: {mode}"
            )
        if len(data) != 112:
            raise ValueError("Expected a 112-byte mode-10 TOF calibration")
        if struct.unpack_from("<3i", data, 60) != (3, 68, 2):
            raise NotImplementedError(
                "Unsupported mode-10 TOF calibration layout"
            )
        values = struct.unpack_from("<7d", data, 4)
        if not all(isfinite(value) for value in values):
            raise ValueError("TOF calibration values must be finite")
        if values[1] <= 0 or values[3] <= 0:
            raise ValueError("Invalid TOF time step or coefficient")
        if values[5] != 0:
            raise NotImplementedError(
                "Nonzero higher-order TOF correction is not supported"
            )
        return cls(*values, raw_blob=data)

    def index_to_mz(self, index: float) -> float:
        if not isfinite(index) or index < 0:
            raise ValueError("TOF index must be finite and nonnegative")
        if (
            self.time_step <= 0
            or self.quadratic_coefficient <= 0
            or not all(
                isfinite(value)
                for value in (
                    self.start_time,
                    self.time_step,
                    self.time_offset,
                    self.quadratic_coefficient,
                    self.mass_time_correction,
                    self.mass_offset,
                )
            )
        ):
            raise ValueError("Invalid TOF calibration values")
        if self.higher_order_correction != 0:
            raise NotImplementedError(
                "Nonzero higher-order TOF correction is not supported"
            )
        time = self.start_time + self.time_step * index - self.time_offset
        linear = 1e6 / sqrt(self.quadratic_coefficient)
        discriminant = (
            linear * linear + 4 * self.mass_time_correction * time
        )
        if time < 0 or not isfinite(discriminant) or discriminant < 0:
            raise ValueError(
                f"TOF calibration is undefined at index {index}"
            )
        # Stable root of time = linear * sqrt(m) + correction * m.
        root_mass = 2 * time / (linear + sqrt(discriminant))
        mz = root_mass * root_mass - self.mass_offset
        if not isfinite(mz) or mz < 0:
            raise ValueError(
                f"TOF calibration is undefined at index {index}"
            )
        return mz


def read_tof_calibration(
    sqlite_path: Path,
    scan_id: int,
    *,
    baf_path: Path,
    index_path: Path,
    calibration_path: Path,
    centroid_offset: int | None = None,
    profile_offset: int | None = None,
) -> TofCalibrationRecord:
    if centroid_offset is None and profile_offset is None:
        raise ValueError(f"Scan {scan_id} has no spectrum blocks")
    if not sqlite_path.is_file():
        raise FileNotFoundError(
            f"TOF reading requires the existing Baf2Sql cache: {sqlite_path}"
        )
    uri = sqlite_path.absolute().as_uri()
    if not uri.startswith("file:///"):
        # SQLite needs an empty URI authority for Windows UNC paths.
        uri = "file:////" + uri.removeprefix("file://")
    with closing(
        sqlite3.connect(uri + "?mode=ro", uri=True)
    ) as connection:
        schema = connection.execute(
            "SELECT Value FROM Properties WHERE Key = 'SchemaType'"
        ).fetchone()
        if schema != ("Baf2Sql",):
            raise ValueError("Expected a Baf2Sql SQLite cache")
        sizes = connection.execute(
            "SELECT BAFSize, IDXSize, XTRSize FROM Info"
        ).fetchall()
        expected_sizes = (
            baf_path.stat().st_size,
            index_path.stat().st_size,
            calibration_path.stat().st_size,
        )
        if sizes != [expected_sizes]:
            raise ValueError(
                "SQLite cache does not match the BAF file sizes"
            )
        row = connection.execute(
            "SELECT s.LineIndexId, s.LineMzId, s.LineIntensityId, "
            "s.ProfileMzId, s.ProfileIntensityId, t.Blob "
            "FROM Spectra AS s LEFT JOIN Transformators AS t "
            "ON t.Id = s.TransformatorId WHERE s.Id = ?",
            (scan_id,),
        ).fetchone()
        if row is None:
            raise ValueError(
                f"Scan {scan_id} is absent from the SQLite cache"
            )
        for offset, identities, kinds, label in (
            (centroid_offset, row[:3], (0x10, 0x16, 0x11), "centroid"),
            (profile_offset, row[3:5], (0x43, 0x42), "profile"),
        ):
            if offset is None:
                continue
            for identity, kind in zip(identities, kinds):
                if (
                    not isinstance(identity, int)
                    or identity >> 56 != kind
                    or identity & ((1 << 56) - 1) != offset
                ):
                    raise ValueError(
                        f"SQLite {label} identity does not match scan {scan_id}"
                    )
        blob = row[5]
        if blob is None:
            raise ValueError(f"Scan {scan_id} has no TOF calibration blob")
        if not isinstance(blob, bytes):
            raise TypeError(
                f"Scan {scan_id} calibration must be a binary blob"
            )
        return TofCalibrationRecord.from_bytes(blob)


_PROFILE_CODEC_HEADER = bytes.fromhex(
    "0000001a0004050607090a0b20010203040506070809000000010000000a"
)
_PROFILE_WIDTHS = (0, 4, 5, 6, 7, 9, 10, 11)
_PROFILE_BASES = (0, 1, 17, 49, 113, 241, 753, 1777)


class _ProfileBits:
    def __init__(self, data: bytes, position: int) -> None:
        self.data = data
        self.position = position

    def read(self, count: int) -> int:
        end = self.position + count
        if end > len(self.data) * 8:
            raise EOFError("Incomplete TOF profile bitstream")
        stop = (end + 7) // 8
        word = int.from_bytes(self.data[self.position // 8 : stop], "big")
        self.position = end
        return (word >> (stop * 8 - end)) & ((1 << count) - 1)


def _decode_tof_profile(payload: bytes) -> list[float]:
    if len(payload) < 38 + 3:
        raise EOFError("Incomplete TOF profile compression header")
    magic, count = struct.unpack_from("<2I", payload)
    if magic != 0xEE77 or payload[8:38] != _PROFILE_CODEC_HEADER:
        raise NotImplementedError("Unsupported TOF profile compression")
    if payload[-3:] != bytes(3):
        raise ValueError("Invalid TOF profile block trailer")
    bits = _ProfileBits(payload[:-3], 38 * 8)
    intensities: list[float] = []
    previous = 0
    while len(intensities) < count:
        selector = 0
        while bits.read(1) == 0:
            selector += 1
            if selector > 9:
                raise ValueError(
                    "Invalid TOF profile variable-length code"
                )
        if selector == 0:
            intensities.append(float(previous))
            continue
        if selector == 9:
            run = bits.read(32)
            if run == 0 or run > count - len(intensities):
                raise ValueError("Invalid TOF profile run length")
            intensities.extend([float(previous)] * run)
            continue
        if selector == 8:
            raw = bits.read(32)
            delta = raw if raw < 0x80000000 else raw - 0x100000000
        else:
            raw = bits.read(_PROFILE_WIDTHS[selector] + 1)
            magnitude = (raw >> 1) + _PROFILE_BASES[selector]
            delta = -magnitude if raw & 1 else magnitude
        # The stored profile is uint32; delta accumulation wraps likewise.
        previous = (previous + delta) & 0xFFFFFFFF
        intensities.append(float(previous))
    if bits.read(16) != 1:
        raise ValueError("Invalid TOF profile end marker")
    padding = (-bits.position) % 32
    if bits.read(padding) != 0 or bits.position != len(bits.data) * 8:
        raise ValueError("Invalid TOF profile alignment padding")
    return intensities


def read_tof_profile_block(baf_path: Path, offset: int) -> list[float]:
    with baf_path.open("rb") as file:
        file.seek(offset)
        header = file.read(44)
        if len(header) != 44:
            raise EOFError("Incomplete TOF profile header")
        block_size, tag, header_size = struct.unpack_from("<3I", header)
        if tag != 0xBFA01001:
            raise ValueError("Expected a TOF profile block")
        if (
            header_size != 44
            or struct.unpack_from("<I", header, 20)[0] != 0x00060102
        ):
            raise NotImplementedError("Unsupported TOF profile encoding")
        payload_size = struct.unpack_from("<I", header, 24)[0]
        if (
            block_size < header_size
            or payload_size != block_size - header_size
        ):
            raise ValueError("Invalid TOF profile block sizes")
        if offset + block_size > baf_path.stat().st_size:
            raise EOFError("Incomplete TOF profile block")
        payload = file.read(payload_size)
        if len(payload) != payload_size:
            raise EOFError("Incomplete TOF profile payload")
    return _decode_tof_profile(payload)


def read_tof_centroid_block(
    baf_path: Path, offset: int
) -> tuple[list[float], list[float]]:
    with baf_path.open("rb") as file:
        file.seek(offset)
        header = file.read(52)
        if len(header) != 52:
            raise EOFError("Incomplete TOF centroid header")
        block_size, tag, header_size = struct.unpack_from("<3I", header)
        if tag != 0xBFA01002:
            raise ValueError("Expected a TOF centroid block")
        if (
            header_size != 52
            or struct.unpack_from("<I", header, 20)[0] != 3
        ):
            raise NotImplementedError("Unsupported TOF centroid encoding")
        if block_size < header_size + 40:
            raise ValueError("TOF centroid block is too short")
        if offset + block_size > baf_path.stat().st_size:
            raise EOFError("Incomplete TOF centroid block")
        payload_header = file.read(40)
        if len(payload_header) != 40:
            raise EOFError("Incomplete TOF centroid array header")
        count, index_offset, intensity_offset = struct.unpack_from(
            "<3I", payload_header
        )
        if count != struct.unpack_from("<I", header, 24)[0]:
            raise ValueError("TOF centroid peak counts do not match")
        if (
            index_offset < 40
            or intensity_offset < index_offset + count * 8
            or intensity_offset + count * 4 > block_size - header_size
        ):
            raise ValueError("Invalid TOF centroid array offsets")
        file.seek(offset + header_size + index_offset)
        indices_data = file.read(count * 8)
        file.seek(offset + header_size + intensity_offset)
        intensities_data = file.read(count * 4)
        if (
            len(indices_data) != count * 8
            or len(intensities_data) != count * 4
        ):
            raise EOFError("Incomplete TOF centroid arrays")
        indices = list(struct.unpack(f"<{count}d", indices_data))
        intensities = list(struct.unpack(f"<{count}f", intensities_data))
        if not all(isfinite(value) and value >= 0 for value in indices):
            raise ValueError(
                "TOF centroid indices must be finite and nonnegative"
            )
        if not all(
            isfinite(value) and value >= 0 for value in intensities
        ):
            raise ValueError(
                "TOF intensities must be finite and nonnegative"
            )
        return indices, intensities
