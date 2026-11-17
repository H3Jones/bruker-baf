import struct
from pathlib import Path

import pytest

from bruker_baf import FticrCalibrationRecord
from bruker_baf.readers import read_baf_calibration


def _calibration_payload(mode: int) -> bytes:
    payload = bytearray(260)
    struct.pack_into("<2i", payload, 0, len(payload), 12)
    struct.pack_into("<i", payload, 16, 60)
    fields_offset = 64
    struct.pack_into("<i", payload, fields_offset + 8, mode)
    struct.pack_into("<3d", payload, fields_offset + 12, 2000, 10, 1000)
    struct.pack_into("<i", payload, fields_offset + 40, 2)
    struct.pack_into("<2d", payload, fields_offset + 48, 0, 100)
    return bytes(payload)


def test_converts_indices_to_mz() -> None:
    calibration = FticrCalibrationRecord(5, 0, 0, 1000, 2, 0, 100)
    assert calibration.index_to_mz(0) == 10
    assert calibration.index_to_mz(1) == 20


@pytest.mark.parametrize("mode", [4, 5])
def test_reads_instrument_specific_record(
    tmp_path: Path, mode: int
) -> None:
    payload = _calibration_payload(mode)
    path = tmp_path / "analysis.baf_xtr"
    path.write_bytes(
        struct.pack("<i4si", 0, b"VECT", len(payload) + 24)
        + payload
        + bytes(12)
        + struct.pack("<i4si", 0, bytes(4), 0)
    )
    calibration = read_baf_calibration(str(path))
    assert type(calibration) is FticrCalibrationRecord
    assert calibration == FticrCalibrationRecord.from_bytes(payload)
    assert calibration.mode == mode
    assert calibration.alpha == (2000 if mode == 5 else 0)
    assert calibration.beta == (10 if mode == 5 else -10)
    denominator = 50 - calibration.beta
    expected = 1000 / denominator + calibration.alpha / 1000
    assert calibration.index_to_mz(1) == pytest.approx(expected)
