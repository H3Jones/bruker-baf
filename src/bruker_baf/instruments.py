import struct
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

_TAG_ACQUISITION = 0xBFA03001
_TAG_STRINGS = 0xBFA03003


class InstrumentType(str, Enum):
    TIMS_TOF = "timsTOF"
    FTMS = "FTMS"


class UnsupportedInstrumentError(ValueError):
    pass


@dataclass(frozen=True)
class InstrumentInfo:
    instrument_type: InstrumentType
    model_name: str
    serial_number: str
    part_number: str


def _read_acquisition_strings(
    baf_path: str | Path, header_offset: int
) -> list[bytes]:
    with open(baf_path, "rb") as file:
        file.seek(header_offset)
        root_header = file.read(8)
        if len(root_header) != 8:
            raise EOFError("BAF file ends before the acquisition header")
        record_size, record_tag = struct.unpack("<II", root_header)
        if record_tag != _TAG_ACQUISITION:
            raise ValueError(
                f"Expected acquisition record at offset {header_offset}, "
                f"found tag 0x{record_tag:08x}"
            )
        if record_size < 8:
            raise ValueError("Invalid acquisition record size")
        file.seek(header_offset)
        record = file.read(record_size)
        if len(record) != record_size:
            raise EOFError("BAF acquisition record is truncated")

    tag_bytes = struct.pack("<I", _TAG_STRINGS)
    search_offset = 8
    while True:
        tag_offset = record.find(tag_bytes, search_offset)
        if tag_offset < 4:
            break
        block_offset = tag_offset - 4
        if block_offset + 12 <= len(record):
            block_size, block_tag, header_size = struct.unpack_from(
                "<III", record, block_offset
            )
            block_end = block_offset + block_size
            if (
                block_tag == _TAG_STRINGS
                and header_size >= 12
                and header_size <= block_size
                and (header_size - 12) % 4 == 0
                and block_end <= len(record)
            ):
                block = record[block_offset:block_end]
                offsets = [header_size]
                for position in range(12, header_size, 4):
                    string_offset = struct.unpack_from(
                        "<I", block, position
                    )[0]
                    if string_offset == 0:
                        break
                    if not header_size <= string_offset < block_size:
                        raise ValueError(
                            "Invalid string offset in acquisition header"
                        )
                    offsets.append(string_offset)

                strings = []
                for string_offset in offsets:
                    terminator = block.find(b"\x00", string_offset)
                    if terminator < 0:
                        raise ValueError(
                            "Unterminated string in acquisition header"
                        )
                    strings.append(block[string_offset:terminator])
                return strings
        search_offset = tag_offset + len(tag_bytes)

    raise ValueError(
        "BAF acquisition header does not contain a valid string table"
    )


def read_instrument_info(
    baf_path: str | Path, header_offset: int
) -> InstrumentInfo:
    if header_offset <= 0:
        raise ValueError("Invalid global header offset")

    with open(baf_path, "rb") as file:
        preamble = file.read(header_offset)
    if len(preamble) != header_offset:
        raise EOFError("BAF file ends before the global header section")

    tokens = preamble.split(b"\x00")
    matches = []
    for index, token in enumerate(tokens):
        normalized = token.strip().lower()
        if normalized.startswith(b"timstof"):
            matches.append((index, InstrumentType.TIMS_TOF, token))
        elif normalized.startswith(b"ftms"):
            matches.append((index, InstrumentType.FTMS, token))

    if len(matches) > 1:
        names = ", ".join(
            instrument_type.value for _, instrument_type, _ in matches
        )
        raise UnsupportedInstrumentError(
            f"Multiple instrument markers found in {baf_path}: {names}"
        )
    if not matches:
        raise UnsupportedInstrumentError(
            f"Could not identify a supported instrument in {baf_path}; "
            "expected a timsTOF or FTMS marker before the global header"
        )

    index, instrument_type, _ = matches[0]
    if index == 0 or index + 1 >= len(tokens):
        raise ValueError(
            f"Incomplete instrument information in {baf_path}"
        )

    serial_token = tokens[index - 1]
    part_token = tokens[index + 1]
    if instrument_type is InstrumentType.FTMS and b"." in part_token:
        part_prefix, serial_suffix = part_token.rsplit(b".", maxsplit=1)
        if part_prefix != serial_token or not serial_suffix.isdigit():
            raise ValueError(
                f"Invalid FTMS part/serial value after marker in {baf_path}"
            )
        part_number = part_prefix
        serial_number = serial_suffix
    else:
        serial_number = serial_token
        part_number = part_token

    if not serial_number.isdigit():
        raise ValueError(
            f"Invalid instrument serial number in {baf_path}"
        )
    if not part_number or not all(
        byte in b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        b"abcdefghijklmnopqrstuvwxyz0123456789._-"
        for byte in part_number
    ):
        raise ValueError(
            f"Invalid instrument part number after marker in {baf_path}"
        )

    acquisition_strings = _read_acquisition_strings(
        baf_path, header_offset
    )
    if len(acquisition_strings) < 2 or not acquisition_strings[1]:
        raise ValueError(
            "BAF acquisition string table does not contain an instrument "
            "model name"
        )
    try:
        model_name = acquisition_strings[1].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(
            "BAF acquisition instrument model name is not valid UTF-8"
        ) from error

    return InstrumentInfo(
        instrument_type=instrument_type,
        model_name=model_name,
        serial_number=serial_number.decode("ascii"),
        part_number=part_number.decode("ascii"),
    )
