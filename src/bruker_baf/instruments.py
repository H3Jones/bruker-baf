from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .readers import BafAcquisitionHeader, read_baf_acquisition_header


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


def read_instrument_info(
    baf_path: str | Path,
    header_offset: int,
    acquisition_header: BafAcquisitionHeader | None = None,
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

    if acquisition_header is None:
        acquisition_header = read_baf_acquisition_header(
            baf_path, header_offset
        )
    model_name = acquisition_header.instrument_name
    if model_name is None:
        raise ValueError(
            "BAF acquisition string table does not contain an instrument "
            "model name"
        )

    return InstrumentInfo(
        instrument_type=instrument_type,
        model_name=model_name,
        serial_number=serial_number.decode("ascii"),
        part_number=part_number.decode("ascii"),
    )
