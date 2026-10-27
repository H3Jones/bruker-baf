import argparse

from utils import (
    read_baf_file_header,
    read_baf_index,
    read_baf_metadata_block,
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
        "--output",
        type=str,
        help="The output file to save the extracted parameters",
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

    # fetch param headers
    file_header = read_baf_file_header(
        full_path(base_directory, main_file[0]), offsets["header"]
    )

    metadata = read_baf_metadata_block(
        full_path(base_directory, main_file[0]),
        offsets["metadata"],
        file_header,
    )

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
