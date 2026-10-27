import argparse
from dataclasses import fields

from bruker_baf import BafReader


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Dump the acquisition-header strings from a BAF dataset"
    )
    parser.add_argument(
        "dot_d_directory",
        help="The directory containing the BAF data files",
    )
    args = parser.parse_args()

    reader = BafReader.from_directory(args.dot_d_directory)
    header = reader.acquisition_header
    print("Acquisition header strings:")
    for index, string_field in enumerate(fields(header.strings)):
        value = getattr(header.strings, string_field.name)
        print(f"  {index} ({string_field.name}): {value}")


if __name__ == "__main__":
    main()
