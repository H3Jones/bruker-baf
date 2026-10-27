from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from .instruments import (
    InstrumentInfo,
    InstrumentType,
    read_instrument_info,
)
from .readers import (
    BafIndex,
    BafMetadataValue,
    CalibrationRecord,
    OffsetDict,
    read_baf_calibration,
    read_baf_file_header,
    read_baf_index,
    read_baf_maxima_block,
    read_baf_metadata_block,
    read_baf_spectrum_block,
)


@dataclass
class Spectrum:
    intensities: list[float]
    mass_to_charge: list[float]


class BafReader:
    def __init__(
        self,
        baf_path: Path,
        index_path: Path,
        calibration_path: Path,
    ) -> None:
        self.baf_path = baf_path
        self.index_path = index_path
        self.calibration_path = calibration_path

        self.index: BafIndex = read_baf_index(str(self.index_path))
        self.parameters = read_baf_file_header(
            str(self.baf_path), self.index.header_offset
        )
        self.calibration: CalibrationRecord | None = None

    @classmethod
    def from_directory(cls, directory: str | Path) -> "BafReader":
        directory = Path(directory)
        if (
            directory.is_file()
            and directory.suffix.lower() == ".d"
            and directory.name.lower() == directory.parent.name.lower()
        ):
            directory = directory.parent
        if not directory.is_dir():
            raise NotADirectoryError(
                f"{directory} is not a valid directory"
            )
        baf_path = cls._find_unique(directory, ".baf")
        index_path = cls._find_unique(directory, ".baf_idx")
        calibration_path = cls._find_unique(directory, ".baf_xtr")
        return cls(baf_path, index_path, calibration_path)

    @staticmethod
    def _find_unique(directory: Path, suffix: str) -> Path:
        matches = sorted(
            path
            for path in directory.iterdir()
            if path.is_file() and path.name.lower().endswith(suffix)
        )
        if not matches:
            raise FileNotFoundError(
                f"No file with suffix {suffix} found in {directory}"
            )
        if len(matches) > 1:
            raise ValueError(
                f"Expected one file with suffix {suffix} in {directory}, "
                f"found {len(matches)}"
            )
        return matches[0]

    def read_calibration(self) -> CalibrationRecord:
        if self.calibration is None:
            self.calibration = read_baf_calibration(
                str(self.calibration_path)
            )
        return self.calibration

    @cached_property
    def instrument_info(self) -> InstrumentInfo:
        return read_instrument_info(
            self.baf_path, self.index.header_offset
        )

    @property
    def instrument_type(self) -> InstrumentType:
        return self.instrument_info.instrument_type

    def _single_scan_offsets(self) -> OffsetDict:
        if len(self.index.sections_by_flag) != 1:
            raise NotImplementedError(
                "Reading BAF files with multiple scan groups is not "
                "implemented yet"
            )
        return next(iter(self.index.sections_by_flag.values()))

    def read_intensities(self) -> list[float]:
        offsets = self._single_scan_offsets()
        return read_baf_spectrum_block(
            str(self.baf_path), offsets["spectrum"]
        )

    def read_spectrum(self) -> Spectrum:
        intensities = self.read_intensities()
        calibration = self.read_calibration()
        mass_to_charge = [
            calibration.index_to_mz(i) for i in range(len(intensities))
        ]
        return Spectrum(
            intensities=intensities, mass_to_charge=mass_to_charge
        )

    def read_maxima(self) -> tuple[list[float], tuple[float, ...]]:
        offsets = self._single_scan_offsets()
        if "maxima" not in offsets:
            raise ValueError("This BAF file has no maxima block")
        return read_baf_maxima_block(
            str(self.baf_path),
            offsets["maxima"],
            self.read_calibration(),
        )

    def read_metadata(self) -> list[BafMetadataValue]:
        offsets = self._single_scan_offsets()
        if "metadata" not in offsets:
            raise ValueError("This BAF file has no metadata block")
        return read_baf_metadata_block(
            str(self.baf_path),
            offsets["metadata"],
            self.parameters,
        )
