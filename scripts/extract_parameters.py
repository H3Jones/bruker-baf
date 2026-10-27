import argparse

from bruker_baf import BafReader


def main() -> None:

    parser = argparse.ArgumentParser(
        description="Process A Bruker BAF file in a given .d directory"
    )
    parser.add_argument(
        "dot_d_directory",
        type=str,
        help="The directory containing the data files",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="The output file to save the extracted parameters",
    )
    args = parser.parse_args()

    reader = BafReader.from_directory(args.dot_d_directory)
    metadata = reader.read_metadata()

    print(f"Found metadata with {len(metadata)} entries")

    if args.output:
        import json

        with open(args.output, "w") as f:
            json.dump([entry.to_dict() for entry in metadata], f)

        print(f"Metadata saved to {args.output}")

    else:
        print("No output file specified, skipping saving metadata.")


if __name__ == "__main__":
    main()
