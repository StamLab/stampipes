from __future__ import unicode_literals

import os, sys, logging, re
import requests
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

    parser.add_argument("-b", "--base_dir", dest="base_dir",
        help="The base directory of the flowcell.")
    parser.add_argument("-p", "--processing_file", dest="processing_file",
        help="The processing_file to use as a guide.")
    
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
        help="Only print out planned symlinks.")
    
    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

def create_links(lane, read, base_dir, dry_run = False):

    alignment = lane["alignments"][0]
    
    sample_dir = os.path.join( base_dir, "Project_%s" % lane["project"], "Sample_%s" % lane["samplesheet_name"] )
    sample_name = alignment["sample_name"]
    
    print "\nCreating links for %s %s\n" % (sample_name, read)
    
    os.chdir(sample_dir)
    
    lane1_fastq = glob.glob("%s_%s_*.fastq.gz" % (sample_name, read))
    
    replace = re.compile(r"_L001$")
    L2_sample_name = replace.sub('_L002', sample_name )
    
    lane2_fastq = glob.glob("%s_%s_*.fastq.gz" % (L2_sample_name, read))
    
    lane1_filecount = len(lane1_fastq)
    lane2_filecount = len(lane2_fastq)
    
    for lane2_filenum in range(1, lane2_filecount+1):
        effective_filenum = lane1_filecount + lane2_filenum        
        orig_filename = "%s_%s_%03d.fastq.gz" % (L2_sample_name, read, lane2_filenum)
        new_filename = "%s_%s_%03d.fastq.gz" % (sample_name, read, effective_filenum)
        
        print "Linking %s => %s" % (orig_filename, new_filename)
        
        if not dry_run:
            os.symlink(orig_filename, new_filename)
        

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

    base_dir = poptions.base_dir

    p = json.loads(open(poptions.processing_file, 'r').read())

    for lane in p['libraries']:
        create_links(lane, "R1", base_dir, poptions.dry_run)
        create_links(lane, "R2", base_dir, poptions.dry_run)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
    

