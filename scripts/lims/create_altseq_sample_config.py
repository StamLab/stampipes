import argparse
import json
import logging
from collections import defaultdict

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger(__name__)


def parse_json(filename):
    with open(filename) as f:
        return json.loads(f.read())


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=parse_json)
    parser.add_argument("--output", default="sample_config.tsv")
    return parser


def group_data(processing_info) -> dict:
    """
    group_data tries to estimate what library pools each library belongs to
    Returns dict of tuple keys, values are a list of library numbers
    """
    output = defaultdict(list)
    for lib in processing_info["libraries"]:
        lib_number = lib["library"]
        key = (
            lib["barcode1"]["reverse_sequence"],
            lib["lane"],
        )
        output[key].append(lib_number)

    return output


def to_tsv(label, data):
    lines = ["pool_name\tsample_name\tlane\tbarcode_index"]
    for datum in sorted(
        data, key=lambda d: (d["lane"], d["pool_name"], d["sample_name"])
    ):
        lines.append(
            "\t".join(
                [
                    label + "_" + datum["pool_name"],
                    datum["sample_name"],
                    str(datum["lane"]),
                    datum["barcode_index"],
                ]
            )
        )
    return "\n".join(lines) + "\n"


def get_config_info(processing_data, ds_number: int):
    pass


def construct_config_entries(data: dict) -> [dict]:
    # Maps library number -> (pool_name, barcode1)
    pool_lookup_table = {}
    for pool, values in data["library_pools"].items():
        value = (pool, values["barcode1"])
        for lib_str in values["libraries"]:
            lib_num = int(lib_str.replace("LN", ""))  # Discard the 'LN' prefix
            if lib_num in pool_lookup_table:
                raise ValueError(
                    "Libnum in more than one pool, %s and %s"
                    % (pool_lookup_table[lib_num], value)
                )
            pool_lookup_table[lib_num] = pool

    results = []
    for library in data["libraries"]:
        datum = {
            "barcode_index": library["barcode_index"],
            "sample_name": library["samplesheet_name"],
            "pool_name": pool_lookup_table[library["library"]],
            "lane": library["lane"],
        }
        results.append(datum)
    return results


def main():
    poptions = parser_setup().parse_args()
    label = poptions.data["flowcell"]["label"]

    entries = construct_config_entries(poptions.data)

    tsv = to_tsv(label, entries)
    with open(poptions.output, "w") as f:
        f.write(tsv)


if __name__ == "__main__":
    main()
