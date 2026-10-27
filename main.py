import argparse

from utils import (
    read_baf_calibration,
    read_baf_file_header,
    read_baf_index,
    read_baf_maxima_block,
    read_baf_metadata_block,
    read_baf_spectrum_block,
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


def main() -> None:

    # accept a file path
    parser = argparse.ArgumentParser(
        description="Process A Bruker BAF file in a given .d directory"
    )
    parser.add_argument(
        "dot_d_directory",
        type=str,
        help="The directory containing the data files",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Whether to plot the results",
    )
    args = parser.parse_args()
    base_directory = args.dot_d_directory

    # .d is actually a directory containing data files
    files = list_files(base_directory)
    print("Files in directory:", files)

    # main file -> .baf
    # index file -> .baf_idx
    # calibration file -> .baf_xtr

    main_file = [f for f in files if f.endswith(".baf")]
    index_file = [f for f in files if f.endswith(".baf_idx")]
    calibration_file = [f for f in files if f.endswith(".baf_xtr")]

    print("Main file:", main_file)
    print("Index file:", index_file)
    print("Calibration file:", calibration_file)

    offsets = read_baf_index(full_path(base_directory, index_file[0]))

    calibration_record = read_baf_calibration(
        full_path(base_directory, calibration_file[0])
    )

    if calibration_record is None:
        raise ValueError("No calibration record found")

    print("Calibration record:", calibration_record)

    intensities = read_baf_spectrum_block(
        full_path(base_directory, main_file[0]), offsets
    )

    print(f"Found {len(intensities)} intensity values")

    mz_values = [
        calibration_record.index_to_mz(i) for i in range(len(intensities))
    ]

    # maxima import
    peak_mz, peak_intensities = read_baf_maxima_block(
        full_path(base_directory, main_file[0]),
        offsets["maxima"],
        calibration_record,
    )

    if args.plot:
        plot(
            mz_values,
            intensities,
            peak_mz,
            peak_intensities,
            base_directory,
        )


def plot(mz_values, intensities, peak_mz, peak_intensities, filename):
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.lineplot(x=mz_values, y=intensities)
    sns.scatterplot(x=peak_mz, y=peak_intensities, color="red")
    # add labels and title
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(f"Mass Spectrum for {filename}")
    plt.show()


if __name__ == "__main__":
    main()
