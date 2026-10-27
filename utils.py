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


def read_baf_spectrum_block(
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


def read_baf_maxima_block(
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


@dataclass
class BafHeaderParameter:
    row_index: int
    row_type: int
    parameter_name: str | None
    parameter_group: str | None
    display_name: str | None
    values: str | None
    format: str | None
    unit: str | None


@dataclass
class BafMetadataValue(BafHeaderParameter):
    raw_value: float | str

    @classmethod
    def from_header(
        cls, parameter: BafHeaderParameter, raw_value: float | str
    ) -> "BafMetadataValue":
        return cls(
            row_index=parameter.row_index,
            row_type=parameter.row_type,
            parameter_name=parameter.parameter_name,
            parameter_group=parameter.parameter_group,
            display_name=parameter.display_name,
            values=parameter.values,
            format=parameter.format,
            unit=parameter.unit,
            raw_value=raw_value,
        )

    @property
    def string_value(self) -> str:
        return self.format_string()

    def to_dict(self) -> dict:
        return {
            "row_index": self.row_index,
            "row_type": self.row_type,
            "parameter_name": self.parameter_name,
            "parameter_group": self.parameter_group,
            "display_name": self.display_name,
            "values": self.values,
            "format": self.format,
            "unit": self.unit,
            "raw_value": self.raw_value,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "BafMetadataValue":
        return cls(
            row_index=data["row_index"],
            row_type=data["row_type"],
            parameter_name=data["parameter_name"],
            parameter_group=data["parameter_group"],
            display_name=data["display_name"],
            values=data["values"],
            format=data["format"],
            unit=data["unit"],
            raw_value=data["raw_value"],
        )

    def format_string(self) -> str:
        if isinstance(self.raw_value, str):
            return self.raw_value

        if self.values is not None:
            rounded_number = round(self.raw_value)
            formatted_value = None
            for entry in self.values.split(";"):
                code, separator, label = entry.partition(":")
                if not separator:
                    raise ValueError(
                        f"Invalid enumerated values for row {self.row_index}"
                    )
                try:
                    matches = float(code.strip()) == rounded_number
                except ValueError as error:
                    raise ValueError(
                        f"Invalid enum code for row {self.row_index}: "
                        f"{code!r}"
                    ) from error
                if matches:
                    formatted_value = label
                    break
            if formatted_value is None:
                raise ValueError(
                    f"No enum label for value {rounded_number} on row "
                    f"{self.row_index}"
                )
        else:
            if self.format is None:
                raise ValueError(
                    f"Missing numeric format for row {self.row_index}"
                )
            try:
                if self.format == "%d":
                    formatted_value = self.format % round(self.raw_value)
                else:
                    formatted_value = self.format % self.raw_value
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError(
                    f"Cannot format numeric metadata for row "
                    f"{self.row_index}"
                ) from error

        if self.unit is not None:
            formatted_value = f"{formatted_value} {self.unit}"
        return formatted_value


def _decode_baf_text(data: bytes, context: str) -> str:
    data = data.rstrip(b"\x00")
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        try:
            return data.decode("cp1252")
        except UnicodeDecodeError as error:
            raise ValueError(
                f"{context} contains text that is neither UTF-8 "
                "nor Windows-1252"
            ) from error


def read_baf_file_header(
    file_path: str, header_offset: int
) -> list[BafHeaderParameter]:
    if header_offset < 0:
        raise ValueError("Header offset cannot be negative")

    def read_exact(file, size: int) -> bytes:
        data = file.read(size)
        if len(data) != size:
            raise EOFError(f"Expected {size} bytes, got {len(data)}")
        return data

    with open(file_path, "rb") as file:
        file.seek(header_offset)
        (
            block_size,
            _delimiter,
            _subdelimiter,
            _header_size,
            next_offset,
        ) = struct.unpack("<i4s4sii", read_exact(file, 20))
        if block_size < 20:
            raise ValueError(
                f"Invalid file header block size: {block_size}"
            )
        if next_offset < 20:
            raise ValueError(f"Invalid next-block offset: {next_offset}")

        file.seek(header_offset + next_offset)
        metadata_block_start = file.tell()
        metadata_block_size = struct.unpack("<i", read_exact(file, 4))[0]
        if metadata_block_size < 4:
            raise ValueError(
                f"Invalid metadata block size: {metadata_block_size}"
            )
        file.seek(metadata_block_start + metadata_block_size)
        header_table_start = file.tell()

        (
            table_block_size,
            _table_delimiter,
            _table_header_size,
            row_count,
        ) = struct.unpack("<i4sHH", read_exact(file, 12))
        if table_block_size < 12:
            raise ValueError(
                f"Invalid header table block size: {table_block_size}"
            )
        table_end = header_table_start + table_block_size

        parameters = []
        for _ in range(row_count):
            row_size = struct.unpack("<H", read_exact(file, 2))[0]
            if row_size < 20:
                raise ValueError(f"Invalid header row size: {row_size}")
            row = b"\x00\x00" + read_exact(file, row_size - 2)
            row_index, row_type = struct.unpack_from("<Hi", row, 2)
            field_offsets = struct.unpack_from("<6H", row, 8)

            present_offsets = [
                offset for offset in field_offsets if offset != 0
            ]
            if any(
                offset < 20 or offset >= row_size
                for offset in present_offsets
            ):
                raise ValueError("Header field offset is outside its row")
            if present_offsets != sorted(present_offsets):
                raise ValueError("Header field offsets are out of order")

            fields = []
            for field_index, field_offset in enumerate(field_offsets):
                if field_offset == 0:
                    fields.append(None)
                    continue
                next_field_offset = next(
                    (
                        offset
                        for offset in field_offsets[field_index + 1 :]
                        if offset != 0
                    ),
                    row_size,
                )
                field_data = row[field_offset:next_field_offset]
                fields.append(
                    _decode_baf_text(field_data, f"Header row {row_index}")
                )

            parameters.append(
                BafHeaderParameter(
                    row_index=row_index,
                    row_type=row_type,
                    parameter_name=fields[0],
                    parameter_group=fields[1],
                    display_name=fields[2],
                    values=fields[3],
                    format=fields[4],
                    unit=fields[5],
                )
            )

        if file.tell() > table_end:
            raise ValueError(
                "Header rows extend beyond the header table block"
            )

    return parameters


def read_baf_metadata_block(
    file_path: str,
    metadata_offset: int,
    parameters: list[BafHeaderParameter],
) -> list[BafMetadataValue]:
    if metadata_offset < 0:
        raise ValueError("Metadata offset cannot be negative")

    parameters_by_index = {}
    for parameter in parameters:
        if parameter.row_index in parameters_by_index:
            raise ValueError(
                f"Duplicate header parameter row index: "
                f"{parameter.row_index}"
            )
        parameters_by_index[parameter.row_index] = parameter

    def read_exact(file, size: int) -> bytes:
        data = file.read(size)
        if len(data) != size:
            raise EOFError(f"Expected {size} bytes, got {len(data)}")
        return data

    def read_value_block_header(
        file, block_name: str
    ) -> tuple[int, int, int, int]:
        block_start = file.tell()
        block_size, _delimiter, header_size, value_count = struct.unpack(
            "<i4sii", read_exact(file, 16)
        )
        if header_size < 16 or header_size > block_size:
            raise ValueError(f"Invalid {block_name} block header size")
        if value_count < 0:
            raise ValueError(f"Negative {block_name} value count")
        if block_size < header_size:
            raise ValueError(f"Invalid {block_name} block size")
        file.seek(block_start + header_size)
        return block_start, block_size, header_size, value_count

    def get_parameter(row_index: int) -> BafHeaderParameter:
        try:
            return parameters_by_index[row_index]
        except KeyError as error:
            raise ValueError(
                f"Metadata row {row_index} has no header parameter"
            ) from error

    def make_value(
        row_index: int, raw_value: float | str
    ) -> BafMetadataValue:
        parameter = get_parameter(row_index)
        return BafMetadataValue.from_header(parameter, raw_value)

    with open(file_path, "rb") as file:
        file.seek(metadata_offset)
        (
            metadata_block_size,
            _metadata_delimiter,
            metadata_header_size,
        ) = struct.unpack("<i4si", read_exact(file, 12))
        if (
            metadata_header_size < 12
            or metadata_header_size > metadata_block_size
        ):
            raise ValueError("Invalid metadata block sizes")
        file.seek(metadata_offset + metadata_header_size)

        (
            numeric_start,
            numeric_size,
            numeric_header_size,
            numeric_count,
        ) = read_value_block_header(file, "numeric metadata")
        if numeric_size < numeric_header_size + numeric_count * 10:
            raise ValueError("Numeric metadata values exceed their block")
        numeric_values = []
        for _ in range(numeric_count):
            row_index, number = struct.unpack("<Hd", read_exact(file, 10))
            numeric_values.append(make_value(row_index, number))
        if file.tell() > numeric_start + numeric_size:
            raise ValueError("Numeric metadata values exceed their block")
        file.seek(numeric_start + numeric_size)

        (
            string_start,
            string_size,
            _string_header_size,
            string_count,
        ) = read_value_block_header(file, "string metadata")
        string_values = []
        for _ in range(string_count):
            row_size = struct.unpack("<H", read_exact(file, 2))[0]
            if row_size < 8:
                raise ValueError(
                    f"Invalid metadata string row size: {row_size}"
                )
            row = b"\x00\x00" + read_exact(file, row_size - 2)
            row_index = struct.unpack_from("<H", row, 2)[0]
            declared_length = struct.unpack_from("<i", row, 4)[0]
            string_length = declared_length - 1
            if string_length < 0 or string_length > row_size - 8:
                raise ValueError(
                    f"Invalid metadata string length for row {row_index}"
                )
            value = _decode_baf_text(
                row[8 : 8 + string_length],
                f"Metadata string row {row_index}",
            )
            string_values.append(make_value(row_index, value))
        if file.tell() > string_start + string_size:
            raise ValueError("String metadata values exceed their block")

    return sorted(
        [*numeric_values, *string_values],
        key=lambda metadata: metadata.row_index,
    )
