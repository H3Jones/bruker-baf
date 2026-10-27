import struct
import tempfile
import unittest
from pathlib import Path

from bruker_baf.readers import read_baf_index, read_baf_spectrum_block


def _index_record(
    id_a: int,
    id_b: int,
    flag: int,
    offset: int,
) -> bytes:
    record = bytearray(24)
    record[0:2] = bytes((id_a, id_b))
    struct.pack_into("<i", record, 8, flag)
    struct.pack_into("<i", record, 16, offset)
    return bytes(record)


def _index_vector(payload: bytes) -> bytes:
    return (
        struct.pack("<i4si", 0, b"VECT", len(payload) + 24)
        + payload
        + bytes(12)
    )


class IndexParsingTests(unittest.TestCase):
    def test_index_groups_sections_by_flag(self) -> None:
        index_bytes = b"".join(
            _index_vector(record)
            for record in (
                _index_record(1, 48, 0, 100),
                _index_record(2, 32, 1, 200),
                _index_record(1, 16, 1, 300),
                _index_record(2, 32, 2, 400),
                _index_record(1, 16, 2, 500),
            )
        ) + struct.pack(
            "<i4si", 0, bytes(4), 0
        )

        with tempfile.TemporaryDirectory() as directory:
            index_path = Path(directory) / "analysis.baf_idx"
            index_path.write_bytes(index_bytes)

            index = read_baf_index(str(index_path))

        self.assertEqual(index.header_offset, 100)
        self.assertEqual(
            index.sections_by_flag,
            {
                1: {"metadata": 200, "spectrum": 300},
                2: {"metadata": 400, "spectrum": 500},
            },
        )

    def test_spectrum_discards_three_zero_padding_bytes(self) -> None:
        payload = struct.pack("<2f", 1.25, 2.5) + bytes(3)
        block = struct.pack(
            "<i4si", 12 + len(payload), b"SPEC", 12
        ) + payload

        with tempfile.TemporaryDirectory() as directory:
            spectrum_path = Path(directory) / "analysis.baf"
            spectrum_path.write_bytes(block)

            intensities = read_baf_spectrum_block(str(spectrum_path), 0)

        self.assertEqual(intensities, [1.25, 2.5])


if __name__ == "__main__":
    unittest.main()
