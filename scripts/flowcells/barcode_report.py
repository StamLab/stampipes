#!/usr/bin/env python

import os, sys, logging, re, json, argparse, glob

script_options = {
    "quiet": False,
    "debug": False,
    "processing_file": os.path.join(os.getcwd(), "processing.json"),
    "basedir": os.getcwd(),
    "json": False,
    "threshold": 1000000,
}

def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processing_file", dest="processing_file",
            help="The processing file to use as a guide.")

    parser.add_argument("-b", "--basedir", dest="basedir",
            help="The base flowcell directory")
    parser.add_argument("-j", "--json", action="store_true", 
            help="Write output in JSON")
    parser.add_argument("-t", "--threshold", dest="threshold", type=int,
            help="The minimum cluster count")
    parser.add_argument("-l", "--lane", dest="lane", type=int,
            help="Report details for only the specified lane")
    parser.set_defaults( **script_options )

    return parser

def sum_barcodes(input_files):
    totals = {}
    for infile in input_files:
        f = open(infile, 'r')
        for line in f:
            words = line.split()
            if len(words) != 2:
                continue
            count, barcode = words
            barcode = barcode.replace('+', '-')
            if not barcode in totals:
                totals[barcode] = 0
            totals[barcode] += int(count)
        f.close()

    return(totals)
 
def get_input_files_for_lane(data, lane, basedir):
    globs = []
    files = []
    for l in data['libraries']:
        # TODO: Allow to work with NoIndex samples that span a lane
        if l['lane'] == lane:
            name = "Project_%s/Sample_%s/%s_%s_L%03d_R1_???.barcodes.txt" % (
                    l['project'], l['samplesheet_name'], l['samplesheet_name'],
                    l['realbarcode'], l['lane'] )

            globs.append(os.path.join(basedir, name))
    globs.append("Undetermined_indices/Sample_lane%d/lane%d_Undetermined_L%03d_R1_???.barcodes.txt" % (lane, lane, lane))
    for g in globs:
        files += glob.glob(g)
    return files

def apply_mask(mask, barcode_string):
    orig_barcodes = barcode_string.split('-')
    while len(orig_barcodes) < len(mask):
        orig_barcodes.append(u'')
    barcodes = [ orig_barcodes[i][:l] for (i, l) in enumerate(mask) ]
    return barcodes

def parse_bases_mask(mask_string):
    mask = map(int, re.findall( r"""(?: i ( \d* ) )""", mask_string, re.X | re.I))
    return mask

def get_expected_barcodes(data):
    mask = parse_bases_mask(data['alignment_group']['bases_mask'])
    libraries = {}
    for l in data['libraries']:
        if l['barcode_index'] == "NoIndex":
            barcode = None
        else:
            barcode = '-'.join( apply_mask(mask, l['barcode_index']) )
        lane = l['lane']
        l['realbarcode'] = barcode
        if not lane in libraries:
            libraries[lane] = {}
        libraries[lane][barcode] = l

    return libraries

def merge_actual_and_expected(expected, actual):
    merged = {}
    fields = ['samplesheet_name', 'id', 'purpose']
    if None in expected:
        total_barcodes = sum( actual.values())
        merged['NoIndex'] = {'cluster_count': total_barcodes}
        for field in fields:
            merged['NoIndex'][field] = expected[None][field]
    else:
        for barcode in actual.keys():
            merged[barcode] = { 'cluster_count': actual[barcode] }
            if barcode in expected:
                for field in fields:
                    merged[barcode][field] = expected[barcode][field]

    return merged

def print_stats_txt(stats, threshold):
    for lane in stats.keys():
        lane_stats = stats[lane]
        print "====== Lane %d =====" % (lane)
        for barcode in sorted(lane_stats.keys(), key = lambda x: lane_stats[x]['cluster_count'], reverse = True):
            m = lane_stats[barcode]
            if (m['cluster_count'] >= threshold) or ('samplesheet_name' in m):
                print "%s,%d,%s,%s" % (barcode, m['cluster_count'], m.get('samplesheet_name',''), m.get('purpose',''))
        print ""

def print_stats_json(stats, threshold):
    for lane in stats.keys():
        stats[lane] = { k: v for k, v in stats[lane].iteritems()
                if v['cluster_count'] >= threshold or 'samplesheet_name' in v }
    print json.dumps(stats, sort_keys=True, indent=2)

def main(args = sys.argv):
    parser = parser_setup()
    poptions = parser.parse_args()

    # Read in our processing file
    process_json = open(poptions.processing_file)
    processing_data = json.load(process_json)
    process_json.close()
    expected = get_expected_barcodes( processing_data )

    if poptions.lane:
        lanes = [poptions.lane]
    else:
        lanes = sorted(list(set([l['lane'] for l in processing_data['libraries']])))

    # Get actual barcodes and merge with expected
    compiled_stats = {}
    for lane in lanes:
        barcode_files = get_input_files_for_lane(processing_data, lane, poptions.basedir)
        actual_barcodes = sum_barcodes(barcode_files)
        
        compiled_stats[lane] = merge_actual_and_expected(expected[lane], actual_barcodes)

    # Print out
    if poptions.json:
        print_stats_json(compiled_stats, poptions.threshold)
    else:
        print_stats_txt(compiled_stats, poptions.threshold)

if __name__ == "__main__":
    main()
