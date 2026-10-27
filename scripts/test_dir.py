import argparse
import struct
from pathlib import Path

from bruker_baf import BafReader


def find_dot_d_directories(folder: str | Path) -> list[Path]:
    folder = Path(folder)
    if folder.is_file() and folder.suffix.lower() == ".d":
        folder = folder.parent
    if folder.is_dir() and folder.suffix.lower() == ".d":
        return [folder]
    return sorted(
        path for path in folder.rglob("*.d") if path.is_dir()
    )


def test_reading_dot_d_directories(folder: str | Path) -> None:
    dot_d_dirs = find_dot_d_directories(folder)
    print(f"Found {len(dot_d_dirs)} .d directories in {folder}")
    if not dot_d_dirs:
        raise ValueError(f"No .d directories found in {folder}")

    failures = 0
    for dot_d_dir in dot_d_dirs:
        try:
            reader = BafReader.from_directory(dot_d_dir)
        except (
            EOFError,
            KeyError,
            OSError,
            ValueError,
            struct.error,
        ) as error:
            failures += 1
            print(f"Failed to open {dot_d_dir}: {error}")
            continue

        try:
            instrument = reader.instrument_info
            instrument_result = (
                f"{instrument.instrument_type.value} "
                f"({instrument.model_name}); "
                f"serial={instrument.serial_number}; "
                f"part={instrument.part_number}"
            )
        except (
            EOFError,
            KeyError,
            OSError,
            ValueError,
            struct.error,
        ) as error:
            failures += 1
            instrument_result = f"failed: {error}"

        try:
            calibration = reader.read_calibration()
            calibration_result = f"mode={calibration.mode}"
        except (
            EOFError,
            KeyError,
            OSError,
            NotImplementedError,
            ValueError,
            struct.error,
        ) as error:
            failures += 1
            calibration_result = f"failed: {error}"

        try:
            metadata = reader.read_metadata()
            metadata_result = f"{len(metadata)} entries"
        except (
            EOFError,
            KeyError,
            OSError,
            NotImplementedError,
            ValueError,
            struct.error,
        ) as error:
            failures += 1
            metadata_result = f"failed: {error}"

        try:
            intensities = reader.read_intensities()
            spectrum_result = f"{len(intensities)} intensities"
        except (
            EOFError,
            KeyError,
            OSError,
            NotImplementedError,
            ValueError,
            struct.error,
        ) as error:
            failures += 1
            spectrum_result = f"failed: {error}"

        print(
            f"Tested {dot_d_dir}: instrument {instrument_result}; "
            f"calibration {calibration_result}; "
            f"metadata {metadata_result}; spectrum {spectrum_result}"
        )

    if failures:
        raise SystemExit(
            f"{failures} parsing operation(s) failed across "
            f"{len(dot_d_dirs)} .d dataset(s)"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Test reading .d directories with BafReader"
    )
    parser.add_argument("folder", help="Folder containing .d directories")
    args = parser.parse_args()

    test_reading_dot_d_directories(args.folder)
