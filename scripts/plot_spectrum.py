import argparse

from bruker_baf import BafReader


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
    calibration = reader.read_calibration()
    spectrum = reader.read_spectrum()
    peak_mz, peak_intensities = reader.read_maxima()

    print("Calibration record:", calibration)
    print(f"Found {len(spectrum.intensities)} intensity values")

    if args.plot:
        plot(
            spectrum.mass_to_charge,
            spectrum.intensities,
            peak_mz,
            peak_intensities,
            args.dot_d_directory,
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
