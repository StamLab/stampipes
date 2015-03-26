from __future__ import unicode_literals

import os, sys, logging, re
#import requests
import json
import fileinput
import argparse
import glob

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

script_options = {
    "quiet": False,
    "debug": False,
    "base_dir": os.getcwd(),
    "processing_file": os.path.join(os.getcwd(), "processing.json"),
    "dry_run": False,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-i", "--input-dir", dest="input_dir",
        help="The input directory to use.")
    parser.add_argument("-o", "--output-dir", dest="output_dir",
        help="The output directory to use.")
    parser.add_argument("-p", "--processing_file", dest="processing_file",
        help="The processing_file to use as a guide.")
    
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
        help="Only print out planned symlinks.")
    
    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

def create_links(lane, read, input_basedir, output_basedir, dry_run = False, undetermined = False):

    sample_name = lane["alignments"][0]["sample_name"]
    short_name = lane["samplesheet_name"]

    if undetermined:
        output_dir = os.path.join( output_basedir, "Undetermined_indices", "Sample_lane1")
    else:
        output_dir = os.path.join( output_basedir, "Project_%s" % lane["project"], "Sample_%s" % lane["samplesheet_name"] )

    # if nextseq
    if True:
        short_name = re.sub(r"_", '-', short_name)
        input_dir  = input_basedir
        input_wildcard = os.path.join(input_dir, "%s_S*_L00?_%s_???.fastq.gz" % (short_name, read))
    else: # eventually could be highseq rapid run linking... have to make some changes
        input_dir = os.path.join( input_basedir,  "Project_%s" % lane["project"], "Sample_%s" % lane["samplesheet_name"] )
        input_wildcard = os.path.join(input_dir, "%s_%s_???.fastq.gz" % (sample_name, read))

    if not dry_run and not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    # This will fail if we have the same sample listed multiple times in the
    # samplesheet (run with different barcodes). 
    # But I've never seen that happen.
    input_fastq = sorted(glob.glob(input_wildcard))

    for idx, input_file in enumerate(input_fastq, start=1):
        output_name = "%s_%s_%03d.fastq.gz" % (sample_name, read, idx)
        output_file = os.path.join(output_dir, output_name)

        rel_path = os.path.relpath(input_file, output_dir)

        print "Linking %s => %s" % (rel_path, output_file)
        if not dry_run and not os.path.exists(output_file):
            os.symlink(rel_path, output_file)
    
def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)

    input_dir = poptions.input_dir

    p = json.loads(open(poptions.processing_file, 'r').read())

    for lane in p['libraries']:
        create_links(lane, "R1", input_dir, poptions.output_dir, poptions.dry_run)
        create_links(lane, "R2", input_dir, poptions.output_dir, poptions.dry_run)

    undet_lane = {"alignments":[{"sample_name": "lane1_Undetermined_L001"}], "samplesheet_name": "Undetermined" }
    for read in ['R1', 'R2']:
        create_links(undet_lane, read, input_dir, poptions.output_dir, poptions.dry_run, True)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
    

