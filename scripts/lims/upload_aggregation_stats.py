#!/usr/bin/env python3
import argparse
import logging
import sys

from stamlims_api.rest import setup_api

# Change to logging.DEBUG to see all messages
logging.basicConfig(level=logging.WARN)

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger(__name__)


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

    parser.add_argument("--aggregation", dest="aggregation_id", type=int)
    parser.add_argument(
        "-f", "--counts_file", nargs="+", help="Tab delimited file of count\tvalue"
    )
    parser.add_argument("--spot", dest="spot_file")

    parser.set_defaults(quiet=False, debug=False)

    return parser


# aggregation
def upload_stats(api, aggregation, stats):
    data = [
        {
            "object_id": aggregation,
            "content_type": "aggregation",
            "stats": stats,
        }
    ]
    response = api.post_single_result(url_addition="stat/create", json=data)
    if response is None:
        raise Exception("Upload failed")

    log.info(response)


def upload_file(api, aggregation, counts_file):
    count_content = open(counts_file, "r")

    log.info("uploading {}".format(counts_file))
    stats = {}
    for line in count_content:
        values = line.split()

        if len(values) < 2:
            log.warn("skipping {}".format(values))
            continue

        stat_type_name = values[0]
        value = values[1]

        try:
            float(value)
        except ValueError:
            log.warn(
                "skipping stat-type '{}' with non-numeric value '{}'".format(
                    stat_type_name, value
                )
            )
            continue

        if not stat_type_name:
            log.warn("skipping {}".format(stat_type_name))
            continue
        stats[stat_type_name] = value
        log.debug("{} : {}".format(stat_type_name, value))
    count_content.close()
    upload_stats(api, aggregation, stats)


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

    api = setup_api()

    if poptions.aggregation_id is None:
        log.critical("No --aggregation specified")
        sys.exit(2)
    # aggregation = aggregations.get_aggregation(api, poptions.aggregation_id)
    if poptions.counts_file:
        for count_file in poptions.counts_file:
            upload_file(api, poptions.aggregation_id, count_file)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
