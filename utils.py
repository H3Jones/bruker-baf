import struct
from dataclasses import dataclass
from math import isfinite

type OffsetDict = dict[str, int]


# Read the index file (.baf_idx) and return a dictionary of section offsets.
# The keys are section names and the values are absolute offsets in the file.
def read_baf_index(file_path: str) -> OffsetDict:
    raw = get_raw_vectors(file_path)

    records = [
        IndexRecord.from_bytes(r) for r in raw if len(r) >= 2 and r[1] != 0
    ]

    offsets = {}
    unknown_records = 0
    for record in records:
        name = record.section_name

        if name == "unknown":
            unknown_records += 1
            continue

        if name in offsets:
            raise ValueError(f"Duplicate section name found: {name}")

        offsets[name] = record.absolute_offset

    if unknown_records > 0:
        print(f"Skipped {unknown_records} unknown offset records")

    return offsets


@dataclass
class IndexRecord:
    id_a: int
    id_b: int
    flag: int
    offset: int
    offset_page: int

    @classmethod
    def from_bytes(cls, data: bytes) -> "IndexRecord":
        if len(data) < 24:
            raise ValueError(
                f"Index record is too short: {len(data)} bytes"
            )

        return cls(
            id_a=data[0],
            id_b=data[1],
            flag=struct.unpack_from("<i", data, 8)[0],
            offset=struct.unpack_from("<i", data, 16)[0],
            offset_page=struct.unpack_from("<i", data, 20)[0],
        )

    @property
    def section_name(self) -> str:
        match self.id_a, self.id_b:
            case 1, 48:
                return "header"
            case 2, 32:
                return "metadata"
            case 1, 16:
                return "spectrum"
            case 2, 16:
                return "maxima"
            case _:
                return "unknown"

    @property
    def absolute_offset(self) -> int:
        return absolute_offset(self.offset, self.offset_page)


# Compute the absolute offset in the file given the offset and offset page.
def absolute_offset(offset: int, offset_page: int) -> int:
    unsigned_low = offset & 0xFFFFFFFF
    return (offset_page << 32) + unsigned_low


def get_raw_vectors(filename) -> list[bytes]:
    raw_vectors = []

    def read_exact(file, size):
        data = file.read(size)
        if len(data) != size:
            raise EOFError(f"Expected {size} bytes, got {len(data)}")
        return data

    with open(filename, "rb") as file:
        while True:
            # R: readBin(..., "integer", size=4, endian="little")
            _temp = struct.unpack("<i", read_exact(file, 4))[0]

            # R: readBin(..., "raw", 4)
            _delimiter = read_exact(file, 4)

            # R: readBin(..., "integer", size=4, endian="little")
            block_size = struct.unpack("<i", read_exact(file, 4))[0]

            if block_size == 0:
                break
            if block_size < 24:
                raise ValueError(f"Invalid block size: {block_size}")

            # The record has 24 bytes of framing: 12-byte header + 12-byte trailer.
            payload = read_exact(file, block_size - 24)
            _closer = read_exact(file, 12)
            raw_vectors.append(payload)

    return raw_vectors


def read_baf_main_raw(
    baf_file_path: str, offsets: OffsetDict
) -> list[float]:
    with open(baf_file_path, "rb") as file:
        block_start = offsets["spectrum"]
        file.seek(block_start)

        header = file.read(12)
        if len(header) != 12:
            raise EOFError("Incomplete spectrum block header")

        block_size, _delimiter, header_size = struct.unpack(
            "<i4si", header
        )

        if not 12 <= header_size <= block_size:
            raise ValueError("Invalid spectrum block sizes")

        payload_size = block_size - header_size
        if payload_size % 4 != 0:
            raise ValueError(
                "Spectrum payload is not a whole number of floats"
            )

        file.seek(block_start + header_size)

        payload = file.read(payload_size)
        if len(payload) != payload_size:
            raise EOFError("Incomplete spectrum payload")

        count = payload_size // 4
        intensities = struct.unpack(f"<{count}f", payload)
        return list(intensities)


@dataclass
class CalibrationRecord:
    mode: int
    # ML3
    alpha: float
    # ML2
    beta: float
    # ML1
    freq_high: float
    size: int
    freq_low: float
    freq_width: float

    @classmethod
    def calibration_fields_offset(cls, data: bytes) -> int | None:
        if len(data) < 8:
            return None

        declared_size = struct.unpack_from("<i", data, 0)[0]
        if declared_size != len(data):
            return None

        header_size = struct.unpack_from("<i", data, 4)[0]
        if header_size < 0:
            return None

        block_size_offset = 4 + header_size
        if block_size_offset + 4 > len(data):
            return None

        block_size = struct.unpack_from("<i", data, block_size_offset)[0]
        if block_size < header_size:
            return None

        fields_offset = block_size_offset + block_size - header_size
        if fields_offset + 64 > len(data):
            return None

        return fields_offset

    @classmethod
    def from_bytes(cls, data: bytes) -> "CalibrationRecord":
        offset = cls.calibration_fields_offset(data)
        if offset is None:
            raise ValueError("Unrecognized calibration record layout")

        mode = struct.unpack_from("<i", data, offset + 8)[0]
        if mode not in (4, 5):
            raise ValueError(f"Unsupported calibration mode: {mode}")

        alpha = (
            struct.unpack_from("<d", data, offset + 12)[0]
            if mode == 5
            else 0.0
        )
        beta = struct.unpack_from("<d", data, offset + 20)[0]

        if mode == 4:
            beta = -beta

        freq_high = struct.unpack_from("<d", data, offset + 28)[0]
        size = struct.unpack_from("<i", data, offset + 40)[0]
        freq_low = struct.unpack_from("<d", data, offset + 48)[0]
        freq_width = struct.unpack_from("<d", data, offset + 56)[0]
        if not all(
            isfinite(value)
            for value in (alpha, beta, freq_high, freq_low, freq_width)
        ):
            raise ValueError("Calibration values must be finite")
        if size <= 0 or freq_high == 0:
            raise ValueError("Invalid calibration size or frequency")

        return cls(
            mode=mode,
            alpha=alpha,
            beta=beta,
            freq_high=freq_high,
            size=size,
            freq_low=freq_low,
            freq_width=freq_width,
        )

    def index_to_mz(self, index: float) -> float:
        if self.size <= 0 or self.freq_high == 0:
            raise ValueError("Invalid calibration size or frequency")

        denominator = (
            self.freq_width * (self.size - index) / self.size - self.beta
        )
        if denominator == 0:
            raise ValueError(f"Calibration is undefined at index {index}")

        return self.freq_high / denominator + self.alpha / self.freq_high


def read_baf_calibration(file_path: str) -> CalibrationRecord:
    matches = [
        raw
        for raw in get_raw_vectors(file_path)
        if CalibrationRecord.calibration_fields_offset(raw) is not None
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Expected one calibration record, found {len(matches)}"
        )

    return CalibrationRecord.from_bytes(matches[0])


def read_baf_maxima(
    file_path: str,
    block_start: int,
    calibration: CalibrationRecord,
) -> tuple[list[float], tuple[float, ...]]:

    def read_exact(file, size):
        data = file.read(size)
        if len(data) != size:
            raise EOFError(f"Expected {size} bytes, got {len(data)}")
        return data

    with open(file_path, "rb") as file:
        file.seek(block_start)

        block_size, _delimiter, header_size = struct.unpack(
            "<i4si", read_exact(file, 12)
        )

        if not 12 <= header_size <= block_size:
            raise ValueError("Invalid maxima block sizes")

        payload_start = block_start + header_size
        payload_size = block_size - header_size

        if payload_size < 8:
            raise ValueError("Maxima payload is too short")

        file.seek(payload_start)
        num_peaks, peak_offset = struct.unpack("<ii", read_exact(file, 8))

        if num_peaks < 0:
            raise ValueError("Negative peak count")

        if peak_offset < 8:
            raise ValueError("Peak arrays overlap the peak-list header")

        arrays_size = num_peaks * 16
        if peak_offset + arrays_size > payload_size:
            raise ValueError("Peak arrays extend beyond the maxima block")

        file.seek(payload_start + peak_offset)

        peak_indices = struct.unpack(
            f"<{num_peaks}d",
            read_exact(file, num_peaks * 8),
        )
        peak_intensities = struct.unpack(
            f"<{num_peaks}d",
            read_exact(file, num_peaks * 8),
        )

    peak_mz = [calibration.index_to_mz(index) for index in peak_indices]

    return peak_mz, peak_intensities
