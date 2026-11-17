import argparse
from typing import Literal

from bruker_baf import BafReader


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Read and optionally plot a Bruker BAF spectrum"
    )
    parser.add_argument(
        "--scan-id",
        type=int,
        help="Scan ID to read (required for datasets with multiple scans)",
    )
    parser.add_argument(
        "--mode",
        choices=("profile", "centroid"),
        help="Defaults to FTMS profile or TOF centroid",
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
    reader = BafReader.from_directory(args.dot_d_directory)
    instrument = reader.instrument_info
    acquisition_header = reader.acquisition_header

    print(
        f"Detected {instrument.model_name}: serial "
        f"{instrument.serial_number}, part {instrument.part_number}"
    )
    print(f"Acquisition started at: {acquisition_header.started_utc}")
    run(reader, args.dot_d_directory, args.plot, args.scan_id, args.mode)


def run(
    reader: BafReader,
    dataset: str,
    should_plot: bool,
    scan_id: int | None = None,
    mode: Literal["profile", "centroid"] | None = None,
) -> None:
    calibration = reader.read_calibration(scan_id)
    spectrum = reader.read_spectrum(scan_id, mode=mode)
    selected_scan = scan_id if scan_id is not None else reader.scan_ids[0]
    if "maxima" in reader.index.sections_by_flag[selected_scan]:
        peak_mz, peak_intensities = reader.read_maxima(scan_id)
    else:
        peak_mz, peak_intensities = [], ()

    print("Calibration record:", calibration)
    print(f"Found {len(spectrum.intensities)} intensity values")

    if should_plot:
        plot(
            spectrum.mass_to_charge,
            spectrum.intensities,
            peak_mz,
            peak_intensities,
            dataset,
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
