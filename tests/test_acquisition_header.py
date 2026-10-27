import struct
import tempfile
import unittest
from datetime import UTC, datetime
from pathlib import Path

from bruker_baf.readers import (
    BafAcquisitionStrings,
    read_baf_acquisition_header,
)

_TAG_ACQUISITION = 0xBFA03001
_TAG_STRINGS = 0xBFA03003
_TAG_VARIABLE_DEFS = 0xBFA03002
_TAG_CALIBRATION = 0xBFA00006


class AcquisitionHeaderTests(unittest.TestCase):
    def test_string_fields_are_cleaned_and_missing_fields_default_to_none(
        self,
    ) -> None:
        strings = BafAcquisitionStrings.from_table(
            [" source path ", "solariX", None, "  ", "version"]
        )

        self.assertEqual(strings.source_path, " source path ")
        self.assertEqual(strings.instrument_name, "solariX")
        self.assertIsNone(strings.operator)
        self.assertIsNone(strings.acquisition_software)
        self.assertEqual(strings.software_version, "version")
        self.assertIsNone(strings.sample_name)
        self.assertIsNone(strings.sample_description)

    def test_reads_documented_acquisition_header_fields(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            baf_path = Path(directory) / "analysis.baf"
            baf_path.write_bytes(self._acquisition_record())

            header = read_baf_acquisition_header(baf_path, 0)

        self.assertEqual(
            header.started_utc,
            datetime(2024, 2, 15, 9, 8, 7, 321000, tzinfo=UTC),
        )
        self.assertEqual(header.utc_offset_minutes, -60)
        self.assertEqual(header.analysis_name, "sample-run")
        self.assertEqual(
            header.source_path, r"D:\raw\sample.d\analysis.baf"
        )
        self.assertEqual(header.instrument_name, "solariX")
        self.assertEqual(header.operator, "operator")
        self.assertEqual(header.acquisition_software, "software")
        self.assertEqual(header.software_version, "version")
        self.assertEqual(header.method_name, "method")
        self.assertEqual(header.dotted_version, "version.part")
        self.assertEqual(header.role, "role")
        self.assertEqual(header.sample_name, "sample-name")
        self.assertEqual(header.sample_description, "sample-description")
        self.assertEqual(header.strings.sample_name, "sample-name")
        self.assertEqual(
            header.strings.sample_description,
            "sample-description",
        )
        self.assertEqual(header.string(2), "operator")
        self.assertEqual(
            header.variables[12].name,
            "Collision_Energy_Act",
        )
        self.assertEqual(header.variables[12].group, "MS")
        self.assertEqual(header.variables[12].label, "Collision energy")
        self.assertEqual(header.variables[12].format, "%0.2f")
        self.assertEqual(header.variables[12].unit, "eV")
        self.assertEqual(header.variables[12].kind, 10)
        self.assertEqual(header.reference_calibrations, [b"CAL!"])

    def test_rejects_non_acquisition_record(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            baf_path = Path(directory) / "analysis.baf"
            baf_path.write_bytes(struct.pack("<II", 8, 0xBFA00008))

            with self.assertRaisesRegex(
                ValueError, "Expected acquisition"
            ):
                read_baf_acquisition_header(baf_path, 0)

    def test_reads_sample_values_after_empty_string_slots(self) -> None:
        strings = [
            r"D:\raw\sample.d\analysis.baf",
            "solariX",
            None,
            "Compass ftmsControl",
            "Compass ftmsControl Version.Build",
            "synthetic-method",
            "Compass ftmsControl",
            None,
            "anonymized-sample",
            "synthetic description",
        ]
        with tempfile.TemporaryDirectory() as directory:
            baf_path = Path(directory) / "analysis.baf"
            baf_path.write_bytes(self._acquisition_record(strings))

            header = read_baf_acquisition_header(baf_path, 0)

        self.assertEqual(header.sample_name, "anonymized-sample")
        self.assertEqual(
            header.sample_description,
            "synthetic description",
        )
        self.assertIsNone(header.operator)
        self.assertIsNone(header.role)
        self.assertIsNone(header.strings.operator)
        self.assertIsNone(header.strings.role)
        self.assertEqual(header.string_label(8), "sample_name")

    @classmethod
    def _acquisition_record(
        cls, strings: list[str | None] | None = None
    ) -> bytes:
        if strings is None:
            strings = [
                r"D:\raw\sample.d\analysis.baf",
                "solariX",
                "operator",
                "software",
                "version",
                "method",
                "version.part",
                "role",
                "sample-name",
                "sample-description",
            ]
        strings_block = cls._string_table(strings)
        variables = cls._variable_definitions()
        calibration = (
            struct.pack("<II", 16, _TAG_CALIBRATION) + bytes(4) + b"CAL!"
        )
        date_block = (
            struct.pack("<BBHhHBBB", 15, 2, 2024, -60, 321, 7, 8, 9)
            + b"\x00"
        )
        body = (
            bytes(40)
            + date_block
            + b"sample-run\x00"
            + strings_block
            + variables
            + calibration
        )
        return struct.pack("<II", 8 + len(body), _TAG_ACQUISITION) + body

    @staticmethod
    def _string_table(strings: list[str | None]) -> bytes:
        header_size = 12 + 4 * (len(strings) - 1)
        first_string = strings[0]
        if first_string is None:
            raise ValueError("The first string-table entry cannot be empty")
        string_data = bytearray(first_string.encode("utf-8") + b"\x00")
        offsets = []
        for value in strings[1:]:
            if value is None:
                offsets.append(0)
                continue
            offsets.append(header_size + len(string_data))
            string_data.extend(value.encode("utf-8") + b"\x00")
        size = header_size + len(string_data)
        header = struct.pack("<III", size, _TAG_STRINGS, header_size)
        if offsets:
            header += struct.pack(f"<{len(offsets)}I", *offsets)
        return header + string_data

    @staticmethod
    def _variable_definitions() -> bytes:
        fields = [
            "Collision_Energy_Act",
            "MS",
            "Collision energy",
            "",
            "%0.2f",
            "eV",
        ]
        entry_data = bytearray()
        offsets = []
        for value in fields:
            if not value:
                offsets.append(0)
                continue
            offsets.append(20 + len(entry_data))
            entry_data.extend(value.encode("utf-8") + b"\x00")

        entry_size = 20 + len(entry_data)
        entry = (
            struct.pack(
                "<HHBB7H",
                entry_size,
                12,
                10,
                0,
                0,
                *offsets,
            )
            + entry_data
        )
        size = 12 + len(entry)
        return (
            struct.pack(
                "<IIHH",
                size,
                _TAG_VARIABLE_DEFS,
                12,
                1,
            )
            + entry
        )


if __name__ == "__main__":
    unittest.main()
