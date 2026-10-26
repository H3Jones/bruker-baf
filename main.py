# .d is actually a directory containing data files
from utils import (
    read_baf_calibration,
    read_baf_index,
    read_baf_main_raw,
    read_baf_maxima,
)


def list_files(directory) -> list[str]:
    import os

    return [
        f
        for f in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, f))
    ]


def full_path(directory, filename):
    import os

    return os.path.join(directory, filename)


def main():

    # accept a file path
    import sys

    if len(sys.argv) < 2:
        raise ValueError(
            "Please provide the example directory as an argument"
        )
    example_directory = sys.argv[1]

    # .d is actually a directory containing data files
    files = list_files(example_directory)
    print(files)

    # main file -> .baf
    # index file -> .baf_idx
    # calibration file -> .baf_xtr

    main_file = [f for f in files if f.endswith(".baf")]
    index_file = [f for f in files if f.endswith(".baf_idx")]
    calibration_file = [f for f in files if f.endswith(".baf_xtr")]

    print("Main file:", main_file)
    print("Index file:", index_file)
    print("Calibration file:", calibration_file)

    offsets = read_baf_index(full_path(example_directory, index_file[0]))

    calibration_record = read_baf_calibration(
        full_path(example_directory, calibration_file[0])
    )

    if calibration_record is None:
        raise ValueError("No calibration record found")

    print("Calibration record:", calibration_record)

    intensities = read_baf_main_raw(
        full_path(example_directory, main_file[0]), offsets
    )

    print(f"Found {len(intensities)} intensities")

    mz_values = [
        calibration_record.index_to_mz(i) for i in range(len(intensities))
    ]

    # maxima import
    peak_mz, peak_intensities = read_baf_maxima(
        full_path(example_directory, main_file[0]),
        offsets["maxima"],
        calibration_record,
    )

    plot(mz_values, intensities, peak_mz, peak_intensities)


def plot(mz_values, intensities, peak_mz, peak_intensities):
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.lineplot(x=mz_values, y=intensities)
    sns.scatterplot(x=peak_mz, y=peak_intensities, color="red")
    plt.show()


if __name__ == "__main__":
    main()
