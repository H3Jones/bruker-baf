from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Literal

from .instruments import (
    InstrumentInfo,
    InstrumentType,
    UnsupportedInstrumentError,
    read_instrument_info,
)
from .readers import (
    BafAcquisitionHeader,
    BafIndex,
    BafMetadataValue,
    FticrCalibrationRecord,
    OffsetDict,
    read_baf_acquisition_header,
    read_baf_calibration,
    read_baf_file_header,
    read_baf_index,
    read_baf_maxima_block,
    read_baf_metadata_block,
    read_baf_spectrum_block,
)
from .tof import (
    TofCalibrationRecord,
    read_tof_calibration,
    read_tof_centroid_block,
)


@dataclass
class Spectrum:
    intensities: list[float]
    mass_to_charge: list[float]


class BafReader(ABC):
    _default_mode: Literal["profile", "centroid"]
    _instrument_type: InstrumentType

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
        if self.instrument_type is not self._instrument_type:
            raise UnsupportedInstrumentError(
                f"{type(self).__name__} cannot read "
                f"{self.instrument_type.value} data"
            )

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
        index = read_baf_index(str(index_path))
        instrument = read_instrument_info(baf_path, index.header_offset)
        reader_type: type[BafReader]
        match instrument.instrument_type:
            case InstrumentType.FTMS:
                reader_type = FticrBafReader
            case InstrumentType.TIMS_TOF:
                reader_type = TofBafReader
            case instrument_type:
                raise UnsupportedInstrumentError(
                    f"No BAF reader for {instrument_type.value}"
                )
        if cls is not BafReader and cls is not reader_type:
            raise UnsupportedInstrumentError(
                f"{cls.__name__} cannot read "
                f"{instrument.instrument_type.value} data"
            )
        return reader_type(baf_path, index_path, calibration_path)

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

    @abstractmethod
    def read_calibration(
        self, scan_id: int | None = None
    ) -> FticrCalibrationRecord | TofCalibrationRecord: ...

    @cached_property
    def acquisition_header(self) -> BafAcquisitionHeader:
        return read_baf_acquisition_header(
            self.baf_path, self.index.header_offset
        )

    @cached_property
    def instrument_info(self) -> InstrumentInfo:
        return read_instrument_info(
            self.baf_path,
            self.index.header_offset,
            self.acquisition_header,
        )

    @property
    def instrument_type(self) -> InstrumentType:
        return self.instrument_info.instrument_type

    @property
    def scan_ids(self) -> tuple[int, ...]:
        return tuple(sorted(self.index.sections_by_flag))

    def _scan_offsets(self, scan_id: int | None) -> tuple[int, OffsetDict]:
        if scan_id is None:
            if len(self.scan_ids) != 1:
                raise NotImplementedError(
                    "This BAF file has multiple scan groups; "
                    "select a scan_id explicitly"
                )
            scan_id = self.scan_ids[0]
        if scan_id not in self.index.sections_by_flag:
            raise ValueError(f"Unknown BAF scan ID: {scan_id}")
        return scan_id, self.index.sections_by_flag[scan_id]

    def _spectrum_mode(
        self, mode: Literal["profile", "centroid"] | None
    ) -> Literal["profile", "centroid"]:
        if mode is None:
            return self._default_mode
        if mode not in ("profile", "centroid"):
            raise ValueError(f"Unknown spectrum mode: {mode}")
        return mode

    @abstractmethod
    def read_intensities(
        self,
        scan_id: int | None = None,
        *,
        mode: Literal["profile", "centroid"] | None = None,
    ) -> list[float]: ...

    def read_spectrum(
        self,
        scan_id: int | None = None,
        *,
        mode: Literal["profile", "centroid"] | None = None,
    ) -> Spectrum:
        mode = self._spectrum_mode(mode)
        if mode == "centroid":
            mz, intensities = self.read_maxima(scan_id)
            return Spectrum(
                intensities=list(intensities), mass_to_charge=mz
            )
        intensities = self.read_intensities(scan_id, mode=mode)
        calibration = self.read_calibration(scan_id)
        mass_to_charge = [
            calibration.index_to_mz(i) for i in range(len(intensities))
        ]
        return Spectrum(
            intensities=intensities, mass_to_charge=mass_to_charge
        )

    @abstractmethod
    def read_maxima(
        self, scan_id: int | None = None
    ) -> tuple[list[float], tuple[float, ...]]: ...

    def read_metadata(
        self, scan_id: int | None = None
    ) -> list[BafMetadataValue]:
        _, offsets = self._scan_offsets(scan_id)
        if "metadata" not in offsets:
            raise ValueError("This BAF file has no metadata block")
        return read_baf_metadata_block(
            str(self.baf_path),
            offsets["metadata"],
            self.parameters,
        )


class FticrBafReader(BafReader):
    _default_mode: Literal["profile", "centroid"] = "profile"
    _instrument_type = InstrumentType.FTMS

    def __init__(
        self,
        baf_path: Path,
        index_path: Path,
        calibration_path: Path,
    ) -> None:
        super().__init__(baf_path, index_path, calibration_path)
        self._calibration: FticrCalibrationRecord | None = None

    def read_calibration(
        self, scan_id: int | None = None
    ) -> FticrCalibrationRecord:
        if scan_id is not None:
            self._scan_offsets(scan_id)
        if self._calibration is None:
            self._calibration = read_baf_calibration(
                str(self.calibration_path)
            )
        return self._calibration

    def read_intensities(
        self,
        scan_id: int | None = None,
        *,
        mode: Literal["profile", "centroid"] | None = None,
    ) -> list[float]:
        mode = self._spectrum_mode(mode)
        _, offsets = self._scan_offsets(scan_id)
        if mode == "centroid":
            return list(self.read_maxima(scan_id)[1])
        return read_baf_spectrum_block(
            str(self.baf_path), offsets["spectrum"]
        )

    def read_maxima(
        self, scan_id: int | None = None
    ) -> tuple[list[float], tuple[float, ...]]:
        _, offsets = self._scan_offsets(scan_id)
        if "maxima" not in offsets:
            raise ValueError("This BAF file has no maxima block")
        return read_baf_maxima_block(
            str(self.baf_path),
            offsets["maxima"],
            self.read_calibration(scan_id),
        )


class TofBafReader(BafReader):
    _default_mode: Literal["profile", "centroid"] = "centroid"
    _instrument_type = InstrumentType.TIMS_TOF

    def __init__(
        self,
        baf_path: Path,
        index_path: Path,
        calibration_path: Path,
        sqlite_path: Path | None = None,
    ) -> None:
        super().__init__(baf_path, index_path, calibration_path)
        self.sqlite_path = sqlite_path or baf_path.with_suffix(".sqlite")
        self._calibrations: dict[int, TofCalibrationRecord] = {}

    def read_calibration(
        self, scan_id: int | None = None
    ) -> TofCalibrationRecord:
        scan_id, offsets = self._scan_offsets(scan_id)
        if "maxima" not in offsets:
            raise ValueError(f"Scan {scan_id} has no centroid block")
        if scan_id not in self._calibrations:
            self._calibrations[scan_id] = read_tof_calibration(
                self.sqlite_path.resolve(),
                scan_id,
                baf_path=self.baf_path,
                index_path=self.index_path,
                calibration_path=self.calibration_path,
                centroid_offset=offsets["maxima"],
            )
        return self._calibrations[scan_id]

    def read_intensities(
        self,
        scan_id: int | None = None,
        *,
        mode: Literal["profile", "centroid"] | None = None,
    ) -> list[float]:
        mode = self._spectrum_mode(mode)
        _, offsets = self._scan_offsets(scan_id)
        if mode == "profile":
            raise NotImplementedError(
                "Compressed TOF profile spectra are not supported; "
                "use mode='centroid'"
            )
        if "maxima" not in offsets:
            raise ValueError("This BAF scan has no centroid block")
        return read_tof_centroid_block(self.baf_path, offsets["maxima"])[1]

    def read_maxima(
        self, scan_id: int | None = None
    ) -> tuple[list[float], tuple[float, ...]]:
        _, offsets = self._scan_offsets(scan_id)
        if "maxima" not in offsets:
            raise ValueError("This BAF file has no maxima block")
        calibration = self.read_calibration(scan_id)
        indices, intensities = read_tof_centroid_block(
            self.baf_path, offsets["maxima"]
        )
        return (
            [calibration.index_to_mz(index) for index in indices],
            tuple(intensities),
        )
