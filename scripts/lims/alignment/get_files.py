import os, sys, logging, re
import requests
import json
import fileinput
import argparse
import datetime
import hashlib
import string
from zipfile import ZipFile

token = None
headers = None
lane_tags = None
flowcell_lane_cache = dict()
flowcell_contenttype = None
base_api_url = None

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger("upload_data.py")

script_options = {
    "base_api_url": None,
    "basedir": os.getcwd(),
    "quiet": False,
    "debug": False,
    "alignment_id": None,
    "lane_id": None,
    "library": None,
    "sample": None,
    "sublibrary": None,
    "barcode": None,
    "flowcell": None,
    "lane": None,
    "file_purpose": None,
}


def parser_setup():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Don't print info messages to standard out.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="Print all debug messages to standard out.",
    )

    parser.add_argument(
        "-a",
        "--api",
        dest="base_api_url",
        help="The base API url, if not the default live LIMS.",
    )

    parser.add_argument(
        "-t", "--token", dest="token", help="Your authentication token.  Required."
    )

    parser.add_argument("--aggregation_id", dest="aggregation_id", type=int)
    parser.add_argument("--lane_id", dest="lane_id", type=int)
    parser.add_argument("--lane", dest="lane", type=int)
    parser.add_argument("--library", "--library", dest="library")
    parser.add_argument("--sample", "--sample", dest="sample")
    parser.add_argument("--sublibrary", dest="sublibrary")
    parser.add_argument("--flowcell", dest="flowcell")

    parser.add_argument("-p", "--file_purpose", dest="file_purpose")

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


class FileFetch(object):
    def __init__(self, api_url, token):
        self.api_url = api_url
        self.token = token
        self.headers = {"Authorization": "Token %s" % token}

    def api_single_result(self, url_addition=None, url=None):
        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        request = requests.get(url, headers=self.headers)

        if request.ok:
            logging.debug(request.json())
            return request.json()
        else:
            logging.error("Could not get data from %s" % url)
            logging.error(request)
            return None

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

    def api_single_list_result(self, url_addition=None, url=None, field=None):
        """
        Using a list API url that should bring up a single item, retrieve that single item if it exists.
        """

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        fetch_results = requests.get(url, headers=self.headers)

        if fetch_results.ok:
            results = fetch_results.json()
            if results["count"] > 1:
                log.error("More than one matching item for fetch query: %s" % url)
            elif results["count"] == 0:
                log.debug("No matching items for fetch query: %s" % url)
            else:
                result = results["results"][0]
                log.debug("Single result fetched from %s: %s" % (url, str(result)))
                if field:
                    return result[field]
                return result
        else:
            log.error("Could not execute api query: %s" % url)

        return None

    def get_file_purpose(self, slug):
        filepurpose_url = "file_purpose/?slug=%s" % (slug)

        return self.api_single_list_result(filepurpose_url)

    def get_file_type(self, slug):
        filetype_url = "file_type/?slug=%s" % (slug)

        return self.api_single_list_result(filetype_url)

    def retrieve_file(self, alignment_id, file_purpose):
        alignment = self.api_single_result("flowcell_lane_alignment/%d/" % alignment_id)

        if not alignment:
            logging.critical("Cannot find alignment %d" % alignment_id)
            sys.exit(1)

        logging.debug(alignment)
        generic_relation_query = (
            "object_id=%d&object_content_type=%d&purpose__slug=%s"
            % (alignment["id"], alignment["object_content_type"], file_purpose["slug"])
        )

        files = self.api_list_result("file/?%s" % generic_relation_query)
        if len(files) == 1:
            sys.stdout.write(files[0]["path"] + "\n")
            return

        directories = self.api_list_result("directory/?%s" % generic_relation_query)
        if len(directories) == 1:
            sys.stdout.write(directories[0]["path"] + "\n")

        if len(files) > 1:
            logging.critical(
                "%d %s files found for alignment %d"
                % (len(files), file_purpose["slug"], alignment_id)
            )
            sys.exit(1)
        if len(directories) > 1:
            logging.critical(
                "%d %s directories found for alignment %d"
                % (len(directories), file_purpose["slug"], alignment_id)
            )

        if not files and not directories:
            logging.critical(
                "No files or directories found for alignment %d" % alignment_id
            )
            sys.exit(1)

    def find_single_alignment(self, lane):
        alignments = self.api_list_result(
            "flowcell_lane_alignment/?lane=%d" % lane["id"]
        )

        if len(alignments) > 1:
            logging.warn("More than one alignment found, finding default")

        for alignment in alignments:
            if alignment["default_lane_alignment"]:
                return alignment

        return None

    def find_lanes(self, args):
        query = {}
        if args.flowcell:
            logging.debug("Using flowcell: %s" % args.flowcell)
            if args.flowcell.startswith("FC"):
                args.flowcell = args.flowcell[2:]
            if len(args.flowcell) != 5:
                logging.warn(
                    "Flowcell label %s is not five characters long" % args.flowcell
                )
            query["flowcell__label"] = args.flowcell

        if args.lane_id:
            logging.debug("Using lane id %d" % args.lane_id)
            query["id"] = args.lane_id

        if args.lane:
            logging.debug("Using lane %d" % args.lane)
            query["lane"] = args.lane

        if args.library:
            logging.debug("Using library %s" % args.library)
            library_number = args.library.strip(string.letters)
            try:
                library_number = int(library_number)
            except ValueError:
                logging.critical("Could not turn %s into library number" % args.library)
                sys.exit(1)

            query["library__number"] = library_number

        if args.sample:
            logging.debug("Using sample %s" % args.sample)
            sample_number = args.sample.lstrip(string.letters)
            if sample_number[-1] in string.letters:
                query["library__sub_library"] = sample_number[-1].upper()
                sample_number = sample_number[:-1]
            try:
                sample_number = int(sample_number)
            except ValueError:
                logging.critical("Could not turn %s into sample number" % args.sample)
            query["sample__number"] = sample_number

        return self.api_list_result(
            "flowcell_lane/?%s"
            % "&".join(["%s=%s" % (item, value) for item, value in query.items()])
        )

    def retrieve(self, args):
        file_purpose = self.get_file_purpose(args.file_purpose)

        if not file_purpose:
            logging.critical("Cannot find file purpose %s" % args.file_purpose)
            sys.exit(1)

        if args.alignment_id:
            self.retrieve_file(alignment_id, file_purpose)

        lanes = self.find_lanes(args)

        if not lanes:
            logging.critical("Could not find any lanes with the arguments given")
            sys.exit(1)

        if len(lanes) > 1:
            logging.critical("More than one lane found for arguments given:")
            for lane in lanes:
                logging.error("\t%d [ %s ]" % (lane["id"], lane["view_url"]))
            sys.exit(1)

        alignment = self.find_single_alignment(lanes[0])

        if not alignment:
            logging.critical("Couldn't find an alignment for lane %d" % lanes[0]["id"])
            sys.exit(1)

        self.retrieve_file(alignment["id"], file_purpose)


def main(args=sys.argv):
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

    if not poptions.file_purpose:
        logging.critical("Require file purpose")
        sys.exit(1)

    fetch = FileFetch(api_url, token)

    fetch.retrieve(poptions)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
