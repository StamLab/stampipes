import json
import os
import sys
import argparse
import logging

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')

script_files = {
    "bwa": "%s/processes/bwa/process_bwa_paired_trimmed.bash" % STAMPIPES,
    "fastqc": "%s/processes/fastqc_only.bash" % STAMPIPES,
}

script_contents = {}

script_options = {
    "quiet": False,
    "debug": False,
    "process_config": "processing.json",
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_scriptname": "run.bash",   
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="The master script to run all sample scripts.")
    parser.add_argument("-p", "--process-config", dest="process_config",
        help="The process config to work off of.")
    
    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


class ProcessSetUp(object):
    
    def __init__(self, processing_configfile, qsub_scriptname, outfile):
        
        self.processing_configfile = processing_configfile        
        self.qsub_scriptname = qsub_scriptname
        self.outfile = outfile
        
    def create(self):
        self.processing_scripts = dict()
        self.p = json.loads(open(self.processing_configfile, 'r').read())

        for lane in self.p['libraries']:
            self.create_script(lane)

        self.run_scripts()

    def add_script(self, script_file, samplesheet_name, priority):
        
        if not priority in self.processing_scripts:
            self.processing_scripts[priority] = list()
        
        self.processing_scripts[priority].append((samplesheet_name, script_file))

    def run_scripts(self):

        outfile = open(self.qsub_scriptname, 'w')
        
        for priority in sorted(self.processing_scripts.keys(), reverse=True):
            outfile.write("# Priority %s\n" % str(priority))
            for samplesheet_name, script_file in self.processing_scripts[priority]:
                outfile.write("cd %s && " % os.path.dirname(script_file))
                outfile.write("qsub -N proc%s -cwd -V -S /bin/bash %s\n\n" % (samplesheet_name, self.qsub_scriptname))

        outfile.close()

    def create_script(self, lane):

        alignment = lane["alignments"][0]

        if not alignment['aligner']:
            print "# FastQC only %s" % lane['sample']
            base_script = "fastqc"
        else:
            print "# Aligning %s with bwa" % lane['sample']
            base_script = "bwa"
    
        if not base_script in script_contents:
            script_contents[base_script] = open(script_files[base_script], 'r').read()
    
        script_directory = os.path.join(self.p['alignment_group']['directory'], "Project_%s" % lane['project'], "Sample_%s" % lane['samplesheet_name'])
    
        if not os.path.exists(script_directory):
            print "Creating directory %s" % script_directory
            os.makedirs(script_directory)
    
        script_file = os.path.join( script_directory, self.qsub_scriptname )
        print script_file

        outfile = open(script_file, 'w')
        outfile.write("set -e -o pipefail\n")
        outfile.write("export SAMPLE_NAME=%s\n" % alignment['sample_name'])
        outfile.write("export BWAINDEX=%s\n" % alignment['genome_index_location'])
        outfile.write("export GENOME=%s\n" % alignment['genome_index'])
        outfile.write("export ASSAY=%s\n" % lane['assay'])
        outfile.write("export READLENGTH=%s\n" % self.p['flowcell']['read_length'])
        if self.p['flowcell']['paired_end']:
            outfile.write("export PAIRED=True\n")
        outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane['id'])
        outfile.write("export ALIGNMENT_ID=%s\n" % alignment['id'])
        outfile.write("export FLOWCELL=%s\n" % self.p['flowcell']['label'])
        outfile.write("\n")
        outfile.write(script_contents[base_script])
        outfile.close()
      
        self.add_script(script_file, lane['samplesheet_name'], alignment['priority'])
    
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

    process = ProcessSetUp(poptions.process_config, poptions.sample_scriptname,
        poptions.outfile)

    process.create()

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
