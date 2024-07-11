#!/usr/bin/env python3

import argparse
import json
import logging
import sys

script_options = {
    "processing": "processing.json",
}


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--processing",
        dest="processing",
        help="The JSON file to read barcodes from",
    )

    parser.set_defaults(**script_options)
    return parser


# Generates the barcode reporting mask from processing.json
# Returns a list of all barcode lengths represented in the data
def get_barcode_masks(json_data):
    masks = []
    read_length = json_data["flowcell"]["read_length"]
    index_length = json_data["flowcell"]["index_length"]
    barcode_lengths = get_barcode_lengths(json_data)

    # determines the n's in the mask
    def format_difference(x):
        if index_length - int(x) <= 0:
            return ""
        else:
            return "n" + str(index_length - int(x))

    for length in barcode_lengths:
        parts = length.split("-")

        for i, x in enumerate(parts):
            if x == "0":
                parts[i] = "n" + str(index_length - int(x))
            else:
                parts[i] = "i{}{}".format(x, format_difference(x))

        mask = "y{0},{1},{2},y{0}".format(read_length, parts[0], parts[1])
        masks.append(mask)

    return masks


# Determines how many different sizes of barcodes there are and returns a set of them
def get_barcode_lengths(json_data):
    # in case barcode sequence is null
    max_len = json_data["flowcell"]["index_length"]

    def format_length(x):
        if x:
            return min(max_len, len(x["sequence"]))
        else:
            return 0

    # set of each unique length in the data
    lengths = set(
        [
            "{}-{}".format(
                format_length(lib["barcode1"]), format_length(lib["barcode2"])
            )
            for lib in json_data["libraries"]
        ]
    )

    # Make sure only 1 report is run each for single/dual indexed barcodes until reporting is more flexible
    tempbc1, tempbc2 = [], []
    finalList = []

    for n in lengths:
        if n[2] == "0":
            tempbc1.append(n)
        else:
            tempbc2.append(n)
    if tempbc1 != []:
        finalList.append(sorted(tempbc1)[0])
    if tempbc2 != []:
        finalList.append(sorted(tempbc2)[0])

    return finalList


# Detects if there are barcode collisions per lane
# TODO add support for multiple resolutions, i.e. collision w/ barcode 8bp long and mask i6n2,n8 might not collide @ i8,n8
def detect_collisions(json_data):
    num_lanes = max([lib["lane"] for lib in json_data["libraries"]])

    for i in range(num_lanes):
        barcodes = sorted(
            lib["barcode_index"]
            for lib in json_data["libraries"]
            if lib["lane"] == i + 1 and not lib["failed"]
        )
        if len(barcodes) != len(set(barcodes)):
            collision = [
                barcode
                for x, barcode in enumerate(barcodes)
                if barcode == barcodes[x - 1]
            ]
            logging.error(
                "Collision on lane {}. Barcode(s): {}\n".format(i + 1, collision)
            )
            sys.exit(1)
            return True
    return False


# Returns an array of barcode report masks
# Also detects barcode collisions for each mask
def main():
    """This is the main body of the program that by default uses the arguments
    from the command line."""

    # Parse args
    parser = parser_setup()
    poptions = parser.parse_args()

    # Read JSON
    process_json = open(poptions.processing)
    json_data = json.load(process_json)
    process_json.close()

    # detect_collisions(json_data)
    print(" ".join(get_barcode_masks(json_data)))


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
