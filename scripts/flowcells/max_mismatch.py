#!/usr/bin/env python

# Usage: $0
# Primarily intended for use in other scripts, primarly flowcells/setup.sh
# Algorithm design:
#   Mask barcodes according to use_bases_mask
#   Generate all allowable mismatch settings
#   For each setting (most permissive to least), generate all possible barcodes for each lane
#     If there are any collisions, check the next tighter mismatch setting

import os, sys, logging, re, math
import argparse
import itertools
import json

MAX_MISMATCH_LEVEL = 1  # Nextseq can allow 2, Hiseq 2500 allows only 1
POSSIBLE_MISMATCH_LEVELS = range( MAX_MISMATCH_LEVEL, -1, -1 )

script_options = {
    "processing": "processing.json"
}

def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processing", dest="processing",
            help="The JSON file to read barcodes from")

    parser.set_defaults( **script_options )
    return parser

def gen_snps(word, mismatches):
    for d in range(mismatches+1):
        for locs in itertools.combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in itertools.product(*thisWord):
                yield "".join(poss)

def generate_barcodes(barcode_tuple, mismatch_tuple):
    return set(itertools.product(* [
        gen_snps( barcode_tuple[i], mismatch_tuple[i] )
        for i in range(len(barcode_tuple)) ]))

def is_mismatch_level_okay(barcodes, mismatch_tuple):
    # If there's one or fewer barcodes in a lane, any number of mismatches is okay
    if len(barcodes) <= 1:
        return True

    barcode_collection = set()
    # one barcode in a lane is always okay
    if len(barcodes) == 1:
        return True
    for barcode in barcodes:
        new_barcodes = generate_barcodes( barcode, mismatch_tuple )
        if barcode_collection.isdisjoint( new_barcodes ):
            barcode_collection.update(new_barcodes)
        else:
            return False
    return True

def get_max_mismatch_level(lane_set, index_count):
    mismatch_choices = list( itertools.product( POSSIBLE_MISMATCH_LEVELS, repeat=index_count))
    lanes = lane_set.values() # We don't actually care about lane labels (the key)
    for choice in mismatch_choices:
        no_collisions = all( [
            is_mismatch_level_okay(barcodes, choice)
            for barcodes in lanes ] )
        if no_collisions:
            return choice
    return None

# Takes in mask & barcode, returns list of trimmed barcodes
def apply_mask(mask, barcode_string):
    if barcode_string is None:
        barcode_string = u''
    orig_barcodes = barcode_string.split('-')
    while len(orig_barcodes) < len(mask):
        orig_barcodes.append(u'')
    barcodes = [ orig_barcodes[i][:l] for (i, l) in enumerate(mask) ]
    return barcodes

def create_lane_set(libraries, mask):
    lanes = {}
    for library in libraries:
        lane = library['lane']
        barcodes = tuple(apply_mask(mask, library['barcode_index']))
        if lane not in lanes:
            lanes[lane] = set()
        if barcodes in lanes[lane]:
            sys.stderr.write("Collision on lane %d, barcode %s\n" % ( lane, ','.join(barcodes)))
            sys.exit(1)
        lanes[lane].add(barcodes)
    return lanes

# NB: This assumes that all index reads will start with an i, and be followed by one or more digits
#     e.g: i6n will work, but iiiiiin, i2n3i2, or ni6 will not.
def parse_bases_mask(mask_string):
    mask = map(int, re.findall( r"""(?: i ( \d* ) )""", mask_string, re.X | re.I))
    return mask

def main(args = sys.argv):

    parser = parser_setup()
    poptions = parser.parse_args()

    process_json = open(poptions.processing)
    data = json.load(process_json)
    process_json.close()

    mask_string = data['alignment_group']['bases_mask']
    mask = parse_bases_mask(mask_string)

    # If the flowcell has no index, exit.
    if len(mask) == 0:
        sys.stderr.write("No index reads found, setting mismatches = 1\n")
        print "1"
        sys.exit(0)

    lanes = create_lane_set(data['libraries'], mask)

    mismatch_level = get_max_mismatch_level( lanes, len(mask) )

    if not mismatch_level:
        sys.stderr.write("No allowable mismatch levels found, barcode collision?\n")
        sys.exit(1)

    print ",".join(map(str, mismatch_level))

if __name__ == "__main__":
    main()
