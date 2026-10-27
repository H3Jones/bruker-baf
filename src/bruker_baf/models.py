from dataclasses import dataclass
from pathlib import Path

from .readers import (
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

        self.offsets: OffsetDict = read_baf_index(str(self.index_path))
        self.parameters = None
        self.calibration: CalibrationRecord | None = None

    @classmethod
    def from_directory(cls, directory: str | Path) -> "BafReader":
        directory = Path(directory)
        if not directory.is_dir():
            raise NotADirectoryError(
                f"{directory} is not a valid directory"
            )
        baf_path = cls._find_first(directory, ".baf")
        index_path = cls._find_first(directory, ".baf_idx")
        calibration_path = cls._find_first(directory, ".baf_xtr")
        return cls(baf_path, index_path, calibration_path)

    @staticmethod
    def _find_first(directory: Path, suffix: str) -> Path:
        for file in directory.iterdir():
            if file.suffix == suffix:
                return file
        raise FileNotFoundError(
            f"No file with suffix {suffix} found in {directory}"
        )

    def read_calibration(self) -> CalibrationRecord:
        if self.calibration is None:
            self.calibration = read_baf_calibration(
                str(self.calibration_path)
            )
        return self.calibration

    def _read_raw_intensities(self) -> list[float]:
        return read_baf_spectrum_block(
            str(self.baf_path), self.offsets["spectrum"]
        )

    def read_spectrum(self) -> Spectrum:
        intensities = self._read_raw_intensities()
        calibration = self.read_calibration()
        mass_to_charge = [
            calibration.index_to_mz(i) for i in range(len(intensities))
        ]
        return Spectrum(
            intensities=intensities, mass_to_charge=mass_to_charge
        )

    def read_maxima(self) -> tuple[list[float], tuple[float, ...]]:
        if "maxima" not in self.offsets:
            raise ValueError("This BAF file has no maxima block")
        return read_baf_maxima_block(
            str(self.baf_path),
            self.offsets["maxima"],
            self.read_calibration(),
        )

    def read_metadata(self) -> list[BafMetadataValue]:
        return read_baf_metadata_block(
            str(self.baf_path),
            self.offsets["metadata"],
            read_baf_file_header(
                str(self.baf_path), self.offsets["header"]
            ),
        )
