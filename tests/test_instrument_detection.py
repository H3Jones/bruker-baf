import struct
import tempfile
import unittest
from pathlib import Path

from bruker_baf import (
    InstrumentInfo,
    InstrumentType,
    UnsupportedInstrumentError,
)
from bruker_baf.instruments import read_instrument_info


class InstrumentDetectionTests(unittest.TestCase):
    def test_detects_tims_tof_marker_in_preamble(self) -> None:
        self.assertEqual(
            self._read_info(
                b"AIDA6.4.0.1\x0000000017\x00timsTOF Pro\x0088888888\x00",
                "timsTOF Pro",
            ).instrument_type,
            InstrumentType.TIMS_TOF,
        )

    def test_reads_acquisition_model_and_preamble_serial_and_part(self) -> None:
        info = self._read_info(
            b"AIDA6.8.1\x0000000042\x00FTMS2\x0099999999\x00",
            "solariX",
        )

        self.assertEqual(
            info,
            InstrumentInfo(
                instrument_type=InstrumentType.FTMS,
                model_name="solariX",
                serial_number="00000042",
                part_number="99999999",
            ),
        )

    def test_splits_ftms_part_and_serial_from_composite_token(self) -> None:
        info = self._read_info(
            b"AIDA6.8.1\x0099999999\x00FTMS\x0099999999.00000042\x00",
            "solariX",
        )

        self.assertEqual(info.serial_number, "00000042")
        self.assertEqual(info.part_number, "99999999")

    def test_rejects_unknown_instrument_marker(self) -> None:
        with self.assertRaises(UnsupportedInstrumentError):
            self._read_info(b"AIDA7.0\x00unknown model\x00", "unused")

    @staticmethod
    def _read_info(preamble: bytes, model_name: str) -> InstrumentInfo:
        with tempfile.TemporaryDirectory() as directory:
            baf_path = Path(directory) / "analysis.baf"
            source_path = b"C:\\sample.d\\analysis.baf"
            model = model_name.encode("utf-8")
            strings_header_size = 16
            model_offset = strings_header_size + len(source_path) + 1
            strings_payload = (
                struct.pack(
                    "<IIII",
                    16 + len(source_path) + 1 + len(model) + 1,
                    0xBFA03003,
                    strings_header_size,
                    model_offset,
                )
                + source_path
                + b"\x00"
                + model
                + b"\x00"
            )
            acquisition_record = (
                struct.pack(
                    "<II",
                    48 + len(strings_payload),
                    0xBFA03001,
                )
                + bytes(40)
                + strings_payload
            )
            baf_path.write_bytes(preamble + acquisition_record)
            return read_instrument_info(baf_path, len(preamble))


if __name__ == "__main__":
    unittest.main()
