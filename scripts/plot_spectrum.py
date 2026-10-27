import argparse
from dataclasses import fields

from bruker_baf import (
    BafReader,
    InstrumentType,
    UnsupportedInstrumentError,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Read and optionally plot a Bruker BAF spectrum"
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
    run(reader, args.dot_d_directory, args.plot)


def run(
    reader: BafReader,
    dataset: str,
    should_plot: bool,
) -> None:
    match reader.instrument_type:
        case InstrumentType.TIMS_TOF:
            raise NotImplementedError(
                "Spectrum parsing for TIMS-TOF is not implemented yet"
            )
        case InstrumentType.FTMS:
            pass
        case instrument_type:
            raise UnsupportedInstrumentError(
                f"No calibrated spectrum route for {instrument_type.value}"
            )

    calibration = reader.read_calibration()
    spectrum = reader.read_spectrum()
    peak_mz, peak_intensities = reader.read_maxima()

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
