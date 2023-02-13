#!/usr/bin/env python3

import argparse
import csv
import os
import pathlib
# import pprint
import json

def parse_args():
    parser = argparse.ArgumentParser(
        prog="generate_counts_json.py",
        description="Parses CellRanger-style output and produces a JSON file we can upload to LIMS",
    )
    parser.add_argument("cellranger_directory")
    parser.add_argument("barcode_config_file")
    parser.add_argument("pool_name")
    return parser.parse_args()

def parse_tsv(filename):
    """Parses a TSV with header, return list of dicts"""
    with open(filename) as f:
        return [*csv.DictReader(f, delimiter="\t")]

def parse_linewise_stats(filename):
    """Parses a file with a name-value pair on each line, separated by whitespace"""
    d = {}
    with open(filename) as f:
        for line in f.readlines():
            (key, value) = line.strip().split()
            d[key] = value
    return d

def parse_linewise_csv_stats(filename):
    """Parses a file with a name-value pair on each line, separated by comma"""
    d = {}
    with open(filename) as f:
        for line in f.readlines():
            (key, value) = line.strip().split(",")
            d[key] = value
    return d


def parse_barcode_config(filename):
    """Parse the special barcode config file we use"""
    cfg = {}
    with open(filename) as f:
        for line in f.readlines():
            (sample_name, barcode) = line.strip().split("\t")
            (_illumina_barcode, cell_barcode) = barcode.split("-")
            cfg[cell_barcode] = sample_name
    return cfg

def modify_sample_info(info):
    """ Rewrite the sample stats a bit """
    # Keys to delete from the table
    deletes = [
        "CB",
    ]
    # Keys to rename
    renames = [
        ("cbMatch", "total"),
    ]
    out = info.copy()
    for to_delete in deletes:
        del out[to_delete]
    for old, new in renames:
        out[new] = out[old]
        del out[old]
    return out

def get_sample_stats(opts):
    """
    Gets per-sample stats from the CellReads.stats file
    """
    cfg = parse_barcode_config(opts.barcode_config_file)
    cellreads_path = os.path.join(opts.cellranger_directory, "CellReads.stats")
    sample_counts = parse_tsv(cellreads_path)

    sample_stats = {
        #cfg.get(info['CB']): modify_sample_info(info)
        info['CB']: modify_sample_info(info)
        for info in sample_counts
        if info['CB'] in cfg
    }
    #del sample_stats[None]
    return sample_stats

def get_barcode_stats(opts):
    """ Gets the stats about barcode mapping """
    barcode_path = os.path.join(opts.cellranger_directory, "..", "Barcodes.stats")
    return parse_linewise_stats(barcode_path)

def get_summary_stats(opts):
    """ Gets the Summary stats produced by StarSOLO """
    barcode_path = os.path.join(opts.cellranger_directory, "Summary.csv")
    return parse_linewise_csv_stats(barcode_path)

def get_library_pool_info(opts):
    """ Gets the metadata about the library and pool """
    (flowcell, pool) = opts.pool_name.split("_")
    return {"flowcell_label": flowcell, "pool": pool}

def get_barcode_mapping(opts):
    """ Returns the mapping of barcodes to sample names """
    cfg = parse_barcode_config(opts.barcode_config_file)
    return cfg

def get_all_stats(opts):
    """
    Return all the stats and metadata that this script gathers
    Packaged as a single dict
    """
    pool_info = get_library_pool_info(opts)
    return {
        "barcode_mapping": get_barcode_mapping(opts),
        "barcode_stats": get_barcode_stats(opts),
        "summary_stats": get_summary_stats(opts),
        "samples": get_sample_stats(opts),
        "pool": pool_info["pool"],
        "flowcell_label": pool_info["flowcell_label"],
    }


def main():
    """ Run it all and write to stdout """
    opts = parse_args()
    data = get_all_stats(opts)
    print(json.dumps(data))

if __name__ == "__main__":
    main()
