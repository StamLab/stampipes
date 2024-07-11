import argparse
import json
import sys

MAX_BARCODE_LENGTH = 10


def parseArgs():
    parser = argparse.ArgumentParser(description="Split up fastq files by barcode")
    parser.add_argument(
        "--processing",
        dest="processing_file",
        action="store",
        required=True,
        help="processing.json to use (mandatory)",
    )
    parser.add_argument(
        "--barcodes",
        dest="barcodes_file",
        action="store",
        required=True,
        help="barcode output to compare",
    )
    parser.add_argument(
        "--bcmask",
        dest="barcodes_mask",
        action="store",
        required=True,
        help="barcode mask to assess",
    )
    args = parser.parse_args()
    return args


# Generates the barcode reporting mask from processing.json
# Returns a list of all barcode lengths represented in the data
def get_barcode_masks(json_data):
    masks = []
    read_length = json_data["flowcell"]["read_length"]
    barcode_lengths = get_barcode_lengths(json_data)

    # determines the n's in the mask
    def format_difference(x):
        if MAX_BARCODE_LENGTH - int(x) == 0:
            return ""
        else:
            return "n" + str(MAX_BARCODE_LENGTH - int(x))

    for length in barcode_lengths:
        parts = length.split("-")

        for i, x in enumerate(parts):
            if x == "0":
                parts[i] = "n" + str(MAX_BARCODE_LENGTH - int(x))
            else:
                parts[i] = "i{}{}".format(x, format_difference(x))

        mask = "y{0},{1},{2},y{0}".format(read_length, parts[0], parts[1])
        masks.append(mask)

    return masks


# format length
def format_length(x):
    if x:
        return len(x["sequence"])
    else:
        return "0"


# Determines how many different sizes of barcodes there are and returns a set of them
def get_barcode_lengths(json_data):
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


def main(argv):
    args = parseArgs()
    barcodes = json.load(open(args.barcodes_file))
    process = json.load(open(args.processing_file))
    mask = args.barcodes_mask

    success = "TRUE"
    barcode_lengths = get_barcode_lengths(process)
    masks = get_barcode_masks(process)

    for lib in process["libraries"]:
        # only look at libraries with the same mask
        length1 = format_length(lib["barcode1"])
        length2 = format_length(lib["barcode2"])
        lengths = str(length1) + "-" + str(length2)
        index = barcode_lengths.index(lengths)
        if mask == masks[index]:
            # check to see if the barcode is represent in the report
            bcs = ""
            if lib["barcode2"] is not None:
                bc1 = lib["barcode1"]["reverse_sequence"]
                bc2 = lib["barcode2"]["sequence"]
                bcs = bc1 + bc2
            else:
                bc1 = lib["barcode1"]["reverse_sequence"]
                bcs = bc1

            lane = lib["lane"]
            for lane in barcodes["Lanes"]:
                if lane == lane["LaneIndex"]:
                    if bcs in lane["Counts"]:
                        next
                    else:
                        print(lib)
                        success = "FALSE"
        else:
            next

    print(success)


if __name__ == "__main__":
    main(sys.argv)
