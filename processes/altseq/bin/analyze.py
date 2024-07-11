#!/usr/bin/env python3

import argparse
import csv
import os
import pathlib
import pprint


def parse_args():
    parser = argparse.ArgumentParser(
        prog="analyze.py",
        description="Parses CellRanger-style output and summarizes by barcode",
    )
    parser.add_argument("cellreads")
    parser.add_argument("barcode_config_file")
    parser.add_argument("output_directory")
    return parser.parse_args()


def parse_barcode_config(filename):
    cfg = {}
    with open(filename) as f:
        for line in f.readlines():
            (name, barcode) = line.strip().split("\t")
            cfg[barcode] = name
    return cfg


def parse_cellreads(filename):
    with open(filename) as f:
        return [*csv.DictReader(f, delimiter="\t")]


def write_sample(output_directory, sample):
    if not sample.get("name", None):
        # Skip barcodes not in our list
        return
    output = os.path.join(output_directory, ("%s.stats.txt" % sample["name"]))
    output_keys = [
        "cbMatch",
        "cbPerfect",
        "exonic",
        "intronic",
        "mito",
        "genomeU",
        "genomeM",
        "featureU",
        "featureM",
        "nGenesUnique",
        "exonicAS",
        "intronicAS",
    ]
    with open(output, "w") as f:
        for key in output_keys:
            if key in sample:
                f.write("%s\t%s\n" % (key, sample[key]))


def main():
    opts = parse_args()
    cfg = parse_barcode_config(opts.barcode_config_file)
    samples = parse_cellreads(opts.cellreads)
    for sample in samples:
        sample["name"] = cfg.get(sample["CB"], None)

    pathlib.Path(opts.output_directory).mkdir(parents=True, exist_ok=True)
    for sample in samples:
        write_sample(opts.output_directory, sample)


if __name__ == "__main__":
    main()
