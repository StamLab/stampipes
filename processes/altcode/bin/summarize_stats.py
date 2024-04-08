#!/usr/bin/env python

import argparse
import csv
import json
import logging
import math
import os
import pathlib
import pprint
import re

from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        prog="summarize_stats.py",
        description="Parses StarSOLO output and summarizes by barcode",
    )
    parser.add_argument("solo_dir")
    parser.add_argument("pool_info_file")
    parser.add_argument("--subdir", default="GeneFull")
    return parser.parse_args()


def parse_pool_info(filename):
    with open(filename) as f:
        json_data = json.load(f)
        return json_data


def parse_cellreads(filename):
    with open(filename) as f:
        data = []
        for row in csv.DictReader(f, delimiter="\t"):
            for (k, v) in row.items():
                try:
                    row[k] = int(v)
                except:
                    pass
            data.append(row)
        return data

def parse_summary_stats(filename):
    with open(filename) as f:
        data = {}
        for line in f:
            (key, orig_val) = line.strip().split(",")
            val = orig_val
            try:
                val = int(val)
            except ValueError:
                try:
                    val = float(val)
                    if math.isnan(val) or math.isinf(val):
                        # If NaN or Inf, leave as string
                        # because the python json module fucks those up
                        val = orig_val
                except ValueError:
                    pass
            data[key] = val
    return data

def parse_barcode_stats(filename):
    with open(filename) as f:
        data = {}
        for line in f:
            match = re.match(r"\s*(\w+)\s*(\d+)\s*", line)
            (key, val) = (match.group(1), match.group(2))
            data[key] = val
    return data

REVCOM = {"A": "T", "T": "A", "C": "G", "G": "C"}
def revcom(bc):
    if bc is None:
        return None
    return "".join(REVCOM[x] for x in reversed(bc))

def summarize_by_library(pool_info, stats):
    """
    Stats is a list of observed cell barcodes & how many we saw / how well they
    mapped.
    Given the stats and pool_info, sum/group them appropriately by barcode
    Right now it's kind of backwards.
    For example, a library with barcode2 = "TTTAAGCG" will contain all cells
    that end with "_CGCTTAAA" (the reverse complement)
    """
    def build_barcode_to_library_lookup(pool_info, stats):
        barcode_to_library = {}
        for lib in pool_info["libraries"]:
            #bc = revcom(lib["barcode2"])
            # Some old backward-compatibility
            # Newer stuff is at the top of the list
            if "sample_barcode" in lib:
                bc = lib["sample_barcode"]
            elif "barcode2" in lib:
                bc = revcom(lib["sample_barcode"])
            elif "additional_information" in lib and "barcode2" in lib["additional_information"]:
                bc = revcom(lib["additional_information"]["sample_barcode"])
            barcode_to_library[bc] = lib["LN#"]
        return barcode_to_library


    # Stub out keys
    data = {
        "barcode_mapping": {},
        "barcode_stats": {},
        "flowcell_label": {},
        "pool": {},
        "libraries": {},
        "summary_stats": {},
    }

    bc_to_library = build_barcode_to_library_lookup(pool_info, stats)
    libraries = {}
    for cell in stats:
        total_bc = cell["CB"]
        if "_" not in total_bc:
            if total_bc not in ["CBnotInPasslist"]:
                logging.warning("Skipping possible barcode %s", total_bc)
            continue
        (_, _1, bc) = total_bc.split("_")
        if bc not in libraries:
            libraries[bc] = defaultdict(int)
        for (k, v) in cell.items():
            if k == "CB":
                continue
            libraries[bc][k] += int(v)
    # Convert back to strings (ew)
    for bc in libraries:
        for (k, v) in libraries[bc].items():
            libraries[bc][k] = str(v)

    pool_set = set(lib["library_pool"] for lib in pool_info["libraries"])
    assert len(pool_set) == 1, "Should have exactly 1 pool, instead: %s" % pool_set
    data["pool"] = pool_set.pop()
    flowcell_set = set(lib["additional_information"]["flowcell"] for lib in pool_info["libraries"])
    assert len(flowcell_set) == 1, "Pool should have exactly 1 flowcell, instead %s" % flowcell_set
    data["flowcell_label"] = flowcell_set.pop()[2:]

    data["barcode_mapping"] = bc_to_library
    data["libraries"] = libraries

    return data

def summarize_by_sample(pool_info, stats):
    """
    Stats is a list of observed cell barcodes & how many we saw / how well they
    mapped.
    Given the stats and pool_info, sum/group them appropriately by barcode
    Right now it's kind of backwards.
    For example, a library with barcode2 = "TTTAAGCG" will contain all cells
    that end with "_CGCTTAAA" (the reverse complement)
    """
    def build_barcode_to_sample_lookup(pool_info, stats):
        barcode_to_sample = {}
        for lib in pool_info["libraries"]:
            bc = revcom(lib["barcode2"])
            barcode_to_sample[bc] = lib["sample"]
        return barcode_to_sample


    # Stub out keys
    data = {
        "barcode_mapping": {},
        "barcode_stats": {},
        "flowcell_label": {},
        "pool": {},
        "samples": {},
        "summary_stats": {},
    }

    bc_to_sample = build_barcode_to_sample_lookup(pool_info, stats)
    samples = {}
    for cell in stats:
        total_bc = cell["CB"]
        if "_" not in total_bc:
            if total_bc not in ["CBnotInPasslist"]:
                logging.warning("Skipping possible barcode %s", total_bc)
            continue
        (_, _1, bc) = total_bc.split("_")
        if bc not in samples:
            samples[bc] = defaultdict(int)
        for (k, v) in cell.items():
            if k == "CB":
                continue
            samples[bc][k] += int(v)
    # Convert back to strings (ew)
    for bc in samples:
        for (k, v) in samples[bc].items():
            samples[bc][k] = str(v)

    pool_set = set(lib["library_pool"] for lib in pool_info["libraries"])
    assert len(pool_set) == 1, "Should have exactly 1 pool, instead: %s" % pool_set
    data["pool"] = pool_set.pop()
    flowcell_set = set(lib["additional_information"]["flowcell"] for lib in pool_info["libraries"])
    assert len(flowcell_set) == 1, "Pool should have exactly 1 flowcell, instead %s" % flowcell_set
    data["flowcell_label"] = flowcell_set.pop()[2:]

    data["barcode_mapping"] = bc_to_sample
    data["samples"] = samples

    return data

def main():
    opts = parse_args()
    cfg = parse_pool_info(opts.pool_info_file)

    gene_dir = os.path.join(opts.solo_dir, opts.subdir)
    cell_reads_filename = os.path.join(gene_dir, "CellReads.stats")

    samples = parse_cellreads(cell_reads_filename)

    data = summarize_by_library(cfg, samples)
    data["summary_stats"] = parse_summary_stats(os.path.join(gene_dir, "Summary.csv"))
    data["barcode_stats"] = parse_barcode_stats(os.path.join(opts.solo_dir, "Barcodes.stats"))

    print(json.dumps(data))

if __name__ == "__main__":
    main()
