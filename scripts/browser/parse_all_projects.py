#!/usr/bin/env python3
# pulls out top 1000 projects and their ids

from __future__ import unicode_literals

import argparse
import logging
import os
import sys

import requests

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

script_options = {
    "base_api_url": None,
    "token": None,
    "outfile": "project_list.txt",
}


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--api",
        dest="base_api_url",
        help="The base API url, if not the default live LIMS.",
    )
    parser.add_argument(
        "-t", "--token", dest="token", help="Your authentication token.  Required."
    )
    parser.add_argument(
        "-o", "--outfile", dest="outfile", help="The outfile to save to."
    )
    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)
    return parser


def get_projects(api_url, token, outfile):
    info = requests.get(
        "%s/project/?page_size=1000" % (api_url),
        headers={"Authorization": "Token %s" % token},
    )

    if info.ok:
        result = info.json()
        out = open(outfile, "w")
        for proj in result["results"]:
            outstring = "%s\t%s\n" % (proj["id"], proj["slug"])
            out.write(outstring)

    else:
        logging.error("info was not found within API")

    return


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
        # Set up the logging levels
        logging.getLogger("requests").setLevel(logging.WARNING)
        logging.basicConfig(level=logging.INFO, format=log_format)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
    else:
        logging.error("Could not find LIMS API URL.\n")
        sys.exit(1)

    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        logging.error("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    get_projects(api_url, token, poptions.outfile)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
