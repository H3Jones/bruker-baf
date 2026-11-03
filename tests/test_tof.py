import sqlite3
import struct
from contextlib import closing
from dataclasses import replace
from pathlib import Path

import pytest

from bruker_baf import (
    BafReader,
    FticrBafReader,
    FticrCalibrationRecord,
    InstrumentInfo,
    InstrumentType,
    TofBafReader,
    UnsupportedInstrumentError,
)
from bruker_baf.tof import TofCalibrationRecord, read_tof_centroid_block


def _calibration_blob() -> bytes:
    return struct.pack(
        "<i7d3i",
        10,
        17575.4,
        0.2,
        310.96225984722844,
        151079.99386901484,
        -0.0006699995522587805,
        0.0,
        0.03066096024742993,
        3,
        68,
        2,
    ) + bytes(40)


def _centroid_block(
    indices: tuple[float, ...], intensities: tuple[float, ...]
) -> bytes:
    count = len(indices)
    header = bytearray(52)
    index_offset = 40
    intensity_offset = index_offset + 8 * count
    block_size = 52 + intensity_offset + 4 * count
    struct.pack_into("<3I", header, 0, block_size, 0xBFA01002, 52)
    struct.pack_into("<2I", header, 20, 3, count)
    array_header = bytearray(40)
    struct.pack_into(
        "<3I", array_header, 0, count, index_offset, intensity_offset
    )
    return (
        header
        + array_header
        + struct.pack(f"<{count}d", *indices)
        + struct.pack(f"<{count}f", *intensities)
    )


def _index_vector(kind: tuple[int, int], flag: int, offset: int) -> bytes:
    payload = bytearray(24)
    payload[:2] = bytes(kind)
    struct.pack_into("<i", payload, 8, flag)
    struct.pack_into("<i", payload, 16, offset)
    return struct.pack("<i4si", 0, b"VECT", 48) + payload + bytes(12)


@pytest.fixture
def tof_reader(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> TofBafReader:
    baf_path = tmp_path / "analysis.baf"
    first = _centroid_block(
        (58136.591397849465, 66467.32283464567), (246.0, 334.0)
    )
    second = _centroid_block((78503.23979591837,), (518.0,))
    preamble = b"AIDA6.4\x0000000001\x00timsTOF Pro\x0000000002\x00"
    source = b"C:\\synthetic.d\\analysis.baf\x00"
    model = b"synthetic TOF\x00"
    strings = (
        struct.pack(
            "<4I",
            16 + len(source) + len(model),
            0xBFA03003,
            16,
            16 + len(source),
        )
        + source
        + model
    )
    acquisition = (
        struct.pack("<2I", 48 + len(strings), 0xBFA03001)
        + bytes(40)
        + strings
    )
    prefix = preamble + acquisition
    baf_path.write_bytes(prefix + first + second)
    offsets = (len(prefix), len(prefix) + len(first))
    index_path = tmp_path / "analysis.baf_idx"
    index_path.write_bytes(
        _index_vector((1, 48), 0, len(preamble))
        + b"".join(
            _index_vector(kind, scan_id, offset)
            for scan_id, offset in enumerate(offsets, 1)
            for kind in ((1, 16), (2, 16))
        )
        + struct.pack("<i4si", 0, bytes(4), 0)
    )
    calibration_path = tmp_path / "analysis.baf_xtr"
    calibration_path.write_bytes(b"synthetic xtr")
    with (
        closing(
            sqlite3.connect(tmp_path / "analysis.sqlite")
        ) as connection,
        connection,
    ):
        connection.executescript(
            "CREATE TABLE Properties (Key TEXT PRIMARY KEY, Value TEXT);"
            "INSERT INTO Properties VALUES ('SchemaType', 'Baf2Sql');"
            "CREATE TABLE Info (BAFSize INTEGER, IDXSize INTEGER, XTRSize INTEGER);"
            "CREATE TABLE Transformators (Id INTEGER PRIMARY KEY, Blob BLOB);"
            "CREATE TABLE Spectra (Id INTEGER PRIMARY KEY, TransformatorId INTEGER,"
            "LineIndexId INTEGER, LineMzId INTEGER, LineIntensityId INTEGER);"
        )
        connection.execute(
            "INSERT INTO Info VALUES (?, ?, ?)",
            tuple(
                path.stat().st_size
                for path in (baf_path, index_path, calibration_path)
            ),
        )
        for scan_id, offset in enumerate(offsets, 1):
            blob = _calibration_blob()
            if scan_id == 2:
                data = bytearray(blob)
                struct.pack_into("<d", data, 52, 1.03066096024742993)
                blob = bytes(data)
            connection.execute(
                "INSERT INTO Transformators VALUES (?, ?)", (scan_id, blob)
            )
            connection.execute(
                "INSERT INTO Spectra VALUES (?, ?, ?, ?, ?)",
                (
                    scan_id,
                    scan_id,
                    (0x10 << 56) | offset,
                    (0x16 << 56) | offset,
                    (0x11 << 56) | offset,
                ),
            )
    monkeypatch.setattr(
        "bruker_baf.models.read_baf_file_header", lambda *args: []
    )
    reader = BafReader.from_directory(tmp_path)
    assert isinstance(reader, TofBafReader)
    return reader


class TestTofCalibration:
    def test_matches_vendor_reference_points(self) -> None:
        calibration = TofCalibrationRecord.from_bytes(_calibration_blob())
        for index, expected in (
            (0.0, 45.00052185039294),
            (58136.591397849465, 126.0816186968955),
            (66467.32283464567, 141.04650330608976),
            (166534.88759398495, 386.35555088388287),
            (304497.4126984127, 923.0220394119789),
            (618832.0, 3004.9904669248613),
        ):
            assert calibration.index_to_mz(index) == pytest.approx(
                expected, rel=0, abs=1e-9
            )
        assert calibration.mode == 10
        assert calibration.raw_blob == _calibration_blob()

    def test_handles_zero_correction_without_division_by_zero(
        self,
    ) -> None:
        calibration = replace(
            TofCalibrationRecord.from_bytes(_calibration_blob()),
            mass_time_correction=0,
        )
        time = calibration.start_time - calibration.time_offset
        expected = (
            calibration.quadratic_coefficient * (time * 1e-6) ** 2
            - calibration.mass_offset
        )
        assert calibration.index_to_mz(0) == pytest.approx(expected)

    @pytest.mark.parametrize("index", [-1.0, float("nan"), float("inf")])
    def test_rejects_invalid_indices(self, index: float) -> None:
        calibration = TofCalibrationRecord.from_bytes(_calibration_blob())
        with pytest.raises(ValueError, match="index must be"):
            calibration.index_to_mz(index)

    def test_rejects_undefined_calibration(self) -> None:
        calibration = replace(
            TofCalibrationRecord.from_bytes(_calibration_blob()),
            mass_time_correction=-1000,
        )
        with pytest.raises(ValueError, match="undefined"):
            calibration.index_to_mz(0)

    @pytest.mark.parametrize(
        ("offset", "format_", "value", "error"),
        [
            (0, "<i", 12, NotImplementedError),
            (4, "<d", float("nan"), ValueError),
            (12, "<d", 0.0, ValueError),
            (28, "<d", -1.0, ValueError),
            (44, "<d", 0.1, NotImplementedError),
            (60, "<i", 4, NotImplementedError),
        ],
    )
    def test_rejects_unsupported_or_invalid_blobs(
        self,
        offset: int,
        format_: str,
        value: float,
        error: type[Exception],
    ) -> None:
        blob = bytearray(_calibration_blob())
        struct.pack_into(format_, blob, offset, value)
        with pytest.raises(error):
            TofCalibrationRecord.from_bytes(bytes(blob))

    @pytest.mark.parametrize("size", [0, 3, 4, 111, 113])
    def test_rejects_truncated_or_extended_blobs(self, size: int) -> None:
        blob = (_calibration_blob() + bytes(1))[:size]
        with pytest.raises(ValueError):
            TofCalibrationRecord.from_bytes(blob)


class TestTofCentroids:
    def test_reads_double_indices_and_float_intensities(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "analysis.baf"
        path.write_bytes(_centroid_block((1.25, 2.5), (246.0, 334.0)))
        assert read_tof_centroid_block(path, 0) == (
            [1.25, 2.5],
            [246.0, 334.0],
        )

    def test_reads_empty_centroid_block(self, tmp_path: Path) -> None:
        path = tmp_path / "analysis.baf"
        path.write_bytes(_centroid_block((), ()))
        assert read_tof_centroid_block(path, 0) == ([], [])

    @pytest.mark.parametrize(
        ("offset", "value", "error"),
        [
            (4, 0, ValueError),
            (8, 44, NotImplementedError),
            (20, 2, NotImplementedError),
            (24, 100, ValueError),
            (56, 0, ValueError),
            (60, 40, ValueError),
            (60, 1000, ValueError),
        ],
    )
    def test_rejects_invalid_array_headers(
        self,
        tmp_path: Path,
        offset: int,
        value: int,
        error: type[Exception],
    ) -> None:
        data = bytearray(_centroid_block((1.0,), (10.0,)))
        struct.pack_into("<I", data, offset, value)
        path = tmp_path / "analysis.baf"
        path.write_bytes(data)
        with pytest.raises(error):
            read_tof_centroid_block(path, 0)

    @pytest.mark.parametrize("size", [10, 60, 100])
    def test_rejects_truncated_blocks(
        self, tmp_path: Path, size: int
    ) -> None:
        path = tmp_path / "analysis.baf"
        path.write_bytes(_centroid_block((1.0,), (10.0,))[:size])
        with pytest.raises(EOFError):
            read_tof_centroid_block(path, 0)

    @pytest.mark.parametrize(
        ("index", "intensity"), [(float("nan"), 1.0), (1.0, float("inf"))]
    )
    def test_rejects_nonfinite_arrays(
        self, tmp_path: Path, index: float, intensity: float
    ) -> None:
        path = tmp_path / "analysis.baf"
        path.write_bytes(_centroid_block((index,), (intensity,)))
        with pytest.raises(ValueError):
            read_tof_centroid_block(path, 0)


class TestTofReader:
    def test_reads_selected_scan_with_its_own_calibration(
        self, tof_reader: TofBafReader
    ) -> None:
        assert tof_reader.scan_ids == (1, 2)
        first = tof_reader.read_spectrum(1)
        second = tof_reader.read_spectrum(2, mode="centroid")
        assert first.mass_to_charge == pytest.approx(
            [126.0816186968955, 141.04650330608976], abs=1e-9
        )
        assert first.intensities == [246.0, 334.0]
        assert second.mass_to_charge == pytest.approx(
            [163.14859195947193], abs=1e-9
        )
        assert second.intensities == [518.0]
        assert tof_reader.read_intensities(2) == [518.0]
        assert tof_reader.read_maxima(2) == (
            second.mass_to_charge,
            (518.0,),
        )
        assert tof_reader.read_calibration(
            2
        ) is tof_reader.read_calibration(2)

    @pytest.mark.parametrize(
        "method",
        [
            "read_spectrum",
            "read_intensities",
            "read_calibration",
            "read_maxima",
        ],
    )
    def test_requires_explicit_scan_selection(
        self, tof_reader: BafReader, method: str
    ) -> None:
        with pytest.raises(NotImplementedError, match="scan_id"):
            getattr(tof_reader, method)()

    def test_rejects_unknown_scan(self, tof_reader: BafReader) -> None:
        with pytest.raises(ValueError, match="Unknown BAF scan ID"):
            tof_reader.read_spectrum(99)

    def test_rejects_profile_instead_of_reading_compressed_bytes(
        self, tof_reader: BafReader
    ) -> None:
        with pytest.raises(NotImplementedError, match="Compressed TOF"):
            tof_reader.read_spectrum(1, mode="profile")
        with pytest.raises(NotImplementedError, match="Compressed TOF"):
            tof_reader.read_intensities(1, mode="profile")

    def test_missing_cache_is_not_created(
        self, tof_reader: TofBafReader
    ) -> None:
        tof_reader.sqlite_path.unlink()
        with pytest.raises(
            FileNotFoundError, match="existing Baf2Sql cache"
        ):
            tof_reader.read_spectrum(1)
        assert not tof_reader.sqlite_path.exists()

    @pytest.mark.parametrize(
        ("query", "message"),
        [
            ("UPDATE Info SET BAFSize = 0", "file sizes"),
            ("UPDATE Spectra SET LineMzId = 0 WHERE Id = 1", "identity"),
            (
                "DELETE FROM Transformators WHERE Id = 1",
                "no TOF calibration",
            ),
            ("DELETE FROM Spectra WHERE Id = 1", "absent"),
            (
                "UPDATE Properties SET Value = 'other'",
                "Baf2Sql SQLite cache",
            ),
        ],
    )
    def test_rejects_inconsistent_cache(
        self, tof_reader: TofBafReader, query: str, message: str
    ) -> None:
        with (
            closing(sqlite3.connect(tof_reader.sqlite_path)) as connection,
            connection,
        ):
            connection.execute(query)
        with pytest.raises(ValueError, match=message):
            tof_reader.read_spectrum(1)

    def test_cache_is_unchanged_by_reading(
        self, tof_reader: TofBafReader
    ) -> None:
        before = tof_reader.sqlite_path.read_bytes()
        tof_reader.read_spectrum(1)
        assert tof_reader.sqlite_path.read_bytes() == before

    def test_unc_uri_uses_empty_authority(
        self, tof_reader: TofBafReader, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        local_uri = tof_reader.sqlite_path.as_uri() + "?mode=ro"
        connect = sqlite3.connect
        calls = []

        def connect_local(
            database: str, *, uri: bool
        ) -> sqlite3.Connection:
            calls.append(database)
            return connect(local_uri, uri=uri)

        monkeypatch.setattr(
            Path,
            "as_uri",
            lambda self: "file://server/share/analysis.sqlite",
        )
        monkeypatch.setattr(sqlite3, "connect", connect_local)
        assert tof_reader.read_spectrum(1).intensities == [246.0, 334.0]
        assert calls == ["file:////server/share/analysis.sqlite?mode=ro"]

    def test_single_scan_needs_no_selection(
        self, tof_reader: BafReader
    ) -> None:
        tof_reader.index.sections_by_flag.pop(2)
        assert tof_reader.read_spectrum().intensities == [246.0, 334.0]

    def test_ftms_profile_behavior_is_preserved(
        self, tof_reader: BafReader, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        instrument = InstrumentInfo(
            InstrumentType.FTMS, "synthetic FTMS", "00000001", "00000002"
        )
        monkeypatch.setattr(
            "bruker_baf.models.read_instrument_info",
            lambda *args: instrument,
        )
        reader = BafReader.from_directory(tof_reader.baf_path.parent)
        assert isinstance(reader, FticrBafReader)
        assert not hasattr(reader, "_calibrations")
        assert not hasattr(reader, "sqlite_path")
        calibration = FticrCalibrationRecord(5, 0, 0, 1000, 2, 0, 100)
        monkeypatch.setattr(
            "bruker_baf.models.read_baf_calibration",
            lambda *args: calibration,
        )
        monkeypatch.setattr(
            "bruker_baf.models.read_baf_spectrum_block",
            lambda *args: [1.0, 2.0],
        )
        assert reader.read_calibration() is calibration
        with pytest.raises(NotImplementedError, match="scan_id"):
            reader.read_spectrum()
        reader.index.sections_by_flag.pop(2)
        spectrum = reader.read_spectrum()
        assert spectrum.intensities == [1.0, 2.0]
        assert spectrum.mass_to_charge == [10.0, 20.0]
        assert reader.read_calibration() is calibration
        assert reader.read_calibration(1) is calibration

        def read_maxima(
            file_path: str,
            offset: int,
            record: FticrCalibrationRecord,
        ) -> tuple[list[float], tuple[float, ...]]:
            assert record is calibration
            return [15.0], (3.0,)

        monkeypatch.setattr(
            "bruker_baf.models.read_baf_maxima_block", read_maxima
        )
        centroid = reader.read_spectrum(1, mode="centroid")
        assert centroid.mass_to_charge == [15.0]
        assert centroid.intensities == [3.0]
        assert reader.read_intensities(1, mode="centroid") == [3.0]

    def test_factory_returns_instrument_specific_reader(
        self, tof_reader: TofBafReader
    ) -> None:
        reader = TofBafReader.from_directory(tof_reader.baf_path.parent)
        assert isinstance(reader, TofBafReader)
        assert isinstance(reader.read_calibration(1), TofCalibrationRecord)
        assert not hasattr(reader, "_calibration")

    def test_base_reader_is_abstract(
        self, tof_reader: TofBafReader
    ) -> None:
        with pytest.raises(TypeError, match="abstract"):
            BafReader(
                tof_reader.baf_path,
                tof_reader.index_path,
                tof_reader.calibration_path,
            )

    def test_rejects_wrong_reader_family(
        self, tof_reader: TofBafReader
    ) -> None:
        with pytest.raises(
            UnsupportedInstrumentError, match="cannot read"
        ):
            FticrBafReader.from_directory(tof_reader.baf_path.parent)
        with pytest.raises(
            UnsupportedInstrumentError, match="cannot read"
        ):
            FticrBafReader(
                tof_reader.baf_path,
                tof_reader.index_path,
                tof_reader.calibration_path,
            )
