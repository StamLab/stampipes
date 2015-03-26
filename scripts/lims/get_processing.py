from __future__ import unicode_literals

import os, sys, logging, re
import requests
import json
import fileinput
import argparse

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

script_options = {
    "base_api_url": None,
    "token": None,
    "quiet": False,
    "debug": False,
    "alignment_group": None,
    "flowcell": None,
    "outfile": "processing.json",
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-a", "--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.")

    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required.")
    parser.add_argument("-f", "--flowcell", dest="flowcell",
        help="The flowcell we want to get processing info for.")
    parser.add_argument("-g", "--alignment-group", dest="alignment_group", type=int,
        help="A specific aligment group to get processing info for.")        

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="The outfile to save to.")
    
    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser
    
def get_processing_info(api_url, token, id, outfile):
    
    info = requests.get("%s/flowcell_lane_alignment_group/%d/processing_information" % (api_url, id),
        headers={'Authorization': "Token %s" % token})
        
    if info.ok:
        result = info.json()
        print "Writing results to %s" % outfile
        with open(outfile, 'w') as output:
            json.dump(result, output, sort_keys=True, indent=4, separators=(',', ': '))
    else:
        sys.stderr.write("Could not find processing info for alignment group %s" % str(id))
    
    return 
    
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

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
    else:
        sys.stderr.write("Could not find LIMS API URL.\n")
        sys.exit(1)

    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        sys.stderr.write("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    if poptions.flowcell:

        print "Getting alignment groups for %s" % poptions.flowcell
        
        alignment_groups = requests.get("%s/flowcell_lane_alignment_group/?flowcell__label=%s" % (api_url, poptions.flowcell),
            headers={'Authorization': "Token %s" % token})
        
        if not alignment_groups.ok:
            sys.stderr.write("Could not get alignment groups for flowcell")
            print alignment_groups
            sys.exit(1)
        
        results = alignment_groups.json()
        if results["count"] == 0:
           sys.stderr.write("Could not find an alignment group for flowcell %s\n" % poptions.flowcell)
           sys.exit(1)
        if results["count"] > 1:
           sys.stderr.write("More than one alignment group found: \n\n%s") % alignment_groups
           sys.exit(1)

        get_processing_info(api_url, token, results['results'][0]["id"], poptions.outfile)

    if poptions.alignment_group:
        get_processing_info(api_url, token, poptions.alignment_group, poptions.outfile)        

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
