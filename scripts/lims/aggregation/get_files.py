import argparse
import logging
import os
import string
import sys

import requests

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
    "aggregation_id": None,
    "library_number": None,
    "sample_number": None,
    "sublibrary": None,
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
    parser.add_argument("-l", "--library_number", dest="library_number")

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

    def retrieve_file(self, aggregation_id, file_purpose):
        aggregation = self.api_single_result("aggregation/%d" % aggregation_id)

        if not aggregation:
            logging.critical("Cannot find aggregation %d" % aggregation_id)
            sys.exit(1)

        logging.debug(aggregation)

        files = self.api_list_result(
            "file/?object_id=%d&object_content_type=%d&purpose__slug=%s"
            % (
                aggregation["id"],
                aggregation["object_content_type"],
                file_purpose["slug"],
            )
        )

        if len(files) > 1:
            logging.critical(
                "%d %s files found for aggregation %d"
                % (len(files), file_purpose["slug"], aggregation_id)
            )
            sys.exit(1)

        if not files:
            logging.critical(
                "%d %s files found for aggregation %d"
                % (len(files), file_purpose["slug"], aggregation_id)
            )
            sys.exit(1)

        print(files[0]["path"])

    def find_single_aggregation(self, aggregations):
        for aggregation in aggregations:
            if aggregation["default_aggregation"]:
                return aggregation

        return None

    def retrieve_library_file(self, library_number, file_purpose):
        library = self.api_single_list_result("library/?number=%d" % library_number)

        if not library:
            logging.critical("Could not find library %d" % library_number)
            sys.exit(1)

        logging.debug(library)

        aggregations = self.api_list_result("aggregation/?library=%d" % (library["id"]))

        if len(aggregations) > 1:
            use_aggregation = self.find_single_aggregation(aggregations)
            if not use_aggregation:
                logging.critical(
                    "More than one aggregation for library %d and no default found, must specify aggregation id"
                    % (library_number)
                )
                logging.critical(
                    "Options: "
                    + ", ".join([aggregation["id"] for aggregation in aggregations])
                )
                return
            else:
                logging.warn(
                    "More than one aggregation for library %d, using default"
                    % (library_number)
                )
        elif len(aggregations) == 0:
            logging.critical(
                "Cannot find aggregations for library %d" % (library_number)
            )
            return
        elif len(aggregations) == 1:
            use_aggregation = aggregations[0]

        self.retrieve_file(use_aggregation["id"], file_purpose)

    def retrieve(self, aggregation_id, library_number, file_purpose_slug):
        file_purpose = self.get_file_purpose(file_purpose_slug)

        if not file_purpose:
            logging.critical("Cannot find file purpose %s" % file_purpose_slug)
            sys.exit(1)

        if aggregation_id:
            self.retrieve_file(aggregation_id, file_purpose)
        if library_number:
            self.retrieve_library_file(library_number, file_purpose)


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

    if poptions.library_number:
        try:
            library_number = int(poptions.library_number.strip(string.ascii_letters))
        except ValueError:
            logging.critical(
                "Could not get library number from %s" % poptions.library_number
            )
            sys.exit()
    else:
        library_number = None

    fetch = FileFetch(api_url, token)

    fetch.retrieve(poptions.aggregation_id, library_number, poptions.file_purpose)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
