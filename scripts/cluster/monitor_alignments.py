import os, sys, logging, re
import requests
import json
import argparse
import datetime
import subprocess
import xml.dom.minidom

token = None
headers = None
lane_tags = None
flowcell_lane_cache = dict()
flowcell_contenttype = None
base_api_url = None

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger('upload_data.py')

script_options = {
    "base_api_url": None,
    "basedir": os.getcwd(),
    "quiet": False,
    "debug": False,
    "cluster": "sge",
}

ALIGN_REGEX = re.compile(r"ALIGN#(\d+)")

def run_command(args):
    return subprocess.check_output(args, stderr=subprocess.STDOUT).decode("utf-8")

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

    parser.add_argument("-c", "--cluster", dest="cluster",
                        help="The type of cluster: SGE (default) or SLURM")

    return parser

class SGEChecker(object):
    QSTAT_COMMAND = ["qstat", "-xml"]
    QHOST_COMMAND = ["qhost"]

    def parse_jobnames(self):
        qstat_xml = run_command(self.QSTAT_COMMAND)
        dom=xml.dom.minidom.parseString(qstat_xml)

        jobnames = dom.getElementsByTagName('JB_name')

        alignments = set()

        for job in jobnames:
            jobname = job.childNodes[0].data
            log.debug(jobname)
            match = ALIGN_REGEX.search(jobname)
            if match:
                alignments.add(int(match.group(1)))

        log.info("Alignment IDs: %s" % alignments)
        return alignments

    def get_host_info(self):
        qhost = run_command(self.QHOST_COMMAND)
        return qhost



class SLURMChecker(object):

    SQUEUE_COMMAND = ["squeue", "-o", "%j", "--noheader"]
    SNODES_COMMAND = ["scontrol", "show", "nodes"]

    def parse_jobnames(self):
        alignments = set()
        jobnames = run_command(self.SQUEUE_COMMAND).split("\n")
        for jobname in jobnames:
            match = ALIGN_REGEX.search(jobname)
            if match:
                alignments.add(int(match.group(1)))
        log.info("Alignment IDs: %s" % alignments)
        return alignments

    def get_host_info(self):
        snodes = run_command(self.SNODES_COMMAND)
        return snodes

class ClusterMonitor(object):

    checker = None

    def __init__(self, api_url, token, cluster_type=None):
        self.api_url = api_url
        self.token = token
        self.headers = {'Authorization': "Token %s" % token}
        if cluster_type is None:
            cluster_type = "sge"
        cluster_type = cluster_type.lower()
        if cluster_type == "sge":
            self.checker = SGEChecker()
        elif cluster_type == "slurm":
            self.checker = SLURMChecker()
        else:
            log.critical("Invalid cluster type %s", cluster_type)
            sys.exit(1)


    def api_list_result(self, url_addition=None, url=None):

        more = True
        results = []

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        while more:

            logging.debug("Fetching more results for query %s" % url)

            request = requests.get(url, headers=self.headers)

            if not request.ok:
                logging.error(request)
                return None
            more_results = request.json()
            results.extend(more_results["results"])
            if more_results["next"]:
                url = more_results["next"]
            else:
                more = False

        return results

    def run(self):
        currently_running = self.checker.parse_jobnames()
        marked_running = self.lims_currently_processing()

        # What alignments are currently running but not marked yet?
        mark_alignments = currently_running - marked_running
        if mark_alignments:
            log.info("Marking alignments as processing: %s" % str(mark_alignments))
            [self.update_processing_status(align_id, True) for align_id in mark_alignments]

        # What alignments are currently marked but not running?
        finished_alignments = marked_running - currently_running
        if finished_alignments:
            log.info("Alignments no longer processing: %s" % str(finished_alignments))
            [self.update_processing_status(align_id, False) for align_id in finished_alignments]

        self.update_host_info()

    def update_processing_status(self, align_id, processing=True):

        patch_url = "%s/flowcell_lane_alignment/%d/" % (self.api_url, align_id)

        update_result = requests.patch(patch_url, headers = self.headers, data={'currently_processing': processing})

        if not update_result.ok:
            log.critical("Could not update alignment %d: %s" % (align_id, str(update_result)))
            return False

        return True

    def lims_currently_processing(self):

        fetch_results = self.api_list_result("flowcell_lane_alignment/?currently_processing=True")
        lims_process_align_ids = set()

        if fetch_results == None:
            log.critical("Could not get list of currently processing alignments")
            sys.exit(1)

        for result in fetch_results:
            lims_process_align_ids.add(result['id'])
        log.info("Currently marked as processing on LIMS: %s" % str(lims_process_align_ids))
        return lims_process_align_ids


    def update_host_info(self):
        host_usage = self.checker.get_host_info()

        host = run_command("hostname").split(".")[0]
        key = "%s-usage" % host
        url = "%s/key_value/?key=%s" % (self.api_url, key)
        key_value = self.get_single_result("%s/key_value/?key=%s" % (self.api_url, key))
        if not key_value:
            log.error("Cannot find \'%s\' key value" % key)
            return

        update = requests.patch(key_value["url"], data={"value": host_usage}, headers=headers)

        if update.ok:
            log.info(update.results)
        else:
            log.error("Could not update %s usage." % host)
            log.error(update.text)



    def get_single_result(self, fetch_url, field=None):
        """
        Using a list API url that should bring up a single item, retrieve that single item if it exists.
        """

        fetch_results = requests.get(fetch_url, headers = self.headers)

        if fetch_results.ok:
            results = fetch_results.json()
            if results['count'] > 1:
                log.error("More than one matching item for fetch query: %s" % fetch_url)
            elif results['count'] == 0:
                log.debug("No matching items for fetch query: %s" % fetch_url)
            else:
                result = results['results'][0]
                log.debug("Single result fetched from %s: %s" % (fetch_url, str(result)))
                if field:
                    return result[field]
                return result
        else:
            log.error("Could not execute api query: %s" % fetch_url)

        return None

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
        # Set up the default logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)
        # Make this a little less noisy by default
        requests_log = logging.getLogger("requests.packages.urllib3.connectionpool")
        requests_log.setLevel(logging.WARN)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
        log.debug("Using LIMS API endpoint: %s from environment" % api_url)
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
        log.debug("Using LIMS API endpoint: %s from options" % api_url)
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

    monitor = ClusterMonitor(api_url, token, cluster_type=poptions.cluster)

    monitor.run()

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
