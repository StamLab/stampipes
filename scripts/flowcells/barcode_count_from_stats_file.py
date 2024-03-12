#!/usr/bin/env python

import os, sys, re
import json
import argparse
import logging

default_options = {
    "min_count": 1000000,
}

def parser_setup():
    parser = argparse.ArgumentParser()
    # Optional
    parser.add_argument("-c", "--min-count", type=int, dest="min_count",
            help="The minimum number of reads to report")
    # Mandatory
    parser.add_argument("-s", "--stats", dest="stats_file",
            required=True,
            help="The JSON file to read stats from. Generally fastq/Stats/Stats.json")
    parser.add_argument("-b", "--basedir", dest="base_dir",
            required=True,
            help="The base directory, like /net/seq/data/sequencers/DATE_A#####_####_FLOWCELL_LABEL")
    parser.add_argument("-m", "--mask", dest="mask",
            required=True,
            help="The barcode mask, like y151,i8,i8,y151")

    parser.set_defaults( **default_options )
    return parser

def main():
    parser = parser_setup()
    poptions = parser.parse_args()
    odata = {
        "Lanes":  [],
        "Mask": poptions.mask,
        "Sequencer": "NovaSeq",
        "BaseDir": poptions.base_dir,
    }
    with open(poptions.stats_file) as f:
        idata = json.load(f)
    
    for lane in idata["UnknownBarcodes"]:
        olane = {
            "LaneIndex": lane["Lane"],
            "Total": None,
            "Pass": None,
            "Counts": {
                bc.replace("+",""): { "Total": count, "Pass": count }
                for (bc, count) in lane["Barcodes"].items()
                if count > poptions.min_count
            },
        }
        olane["Total"] = sum(lane["Barcodes"].values())
        olane["Pass"] = olane["Total"]

        odata["Lanes"].append(olane)

    for conversion_result in idata["ConversionResults"]:
            lane_num = conversion_result["LaneNumber"]
            lane_idx = None
            for (i, olane) in enumerate(odata["Lanes"]):
                if int(olane["LaneIndex"]) == int(lane_num):
                    lane_idx = i
                    break
            if lane_idx is None:
                logging.error("Lane %s not in odata", lane_num)
            for sample_info in conversion_result["DemuxResults"]:
                for metric_info in sample_info["IndexMetrics"]:
                    # Get matching count
                    barcode = metric_info["IndexSequence"].replace("+","")
                    count = metric_info["MismatchCounts"]["0"]
                    # Update out_data
                    odata["Lanes"][lane_idx]["Counts"][barcode] = {"Total": count, "Pass": count}
                    odata["Lanes"][lane_idx]["Total"] += count
                    odata["Lanes"][lane_idx]["Pass"] += count


    print(json.dumps(odata))

if __name__ == "__main__":
    main()
