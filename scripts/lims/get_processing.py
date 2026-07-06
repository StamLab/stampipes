#!/usr/bin/env python3
from __future__ import unicode_literals

import argparse
import json
import logging
import os
import re
import sys

import requests

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

script_options = {
    "base_api_url": None,
    "token": None,
    "quiet": False,
    "debug": False,
    "alignment_group": None,
    "flowcell": None,
    "project": None,
    "outfile": "processing.json",
}


def parser_setup():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Don't logging.info messages to standard out.",
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
    parser.add_argument(
        "-f",
        "--flowcell",
        dest="flowcell",
        help="The flowcell we want to get processing info for.",
    )
    parser.add_argument(
        "-p",
        "--project",
        dest="project",
        help="The project to get processing info for.",
    )
    parser.add_argument(
        "-g",
        "--alignment-group",
        dest="alignment_group",
        type=int,
        help="A specific aligment group to get processing info for.",
    )
    parser.add_argument(
        "-e",
        "--experiment",
        dest="experiment",
        help="The experiment to get processing info for.",
    )

    parser.add_argument(
        "-o", "--outfile", dest="outfile", help="The outfile to save to."
    )

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


def describe_response(response):
    """Summarize an API response for error messages."""
    body = response.text
    content_type = response.headers.get("content-type", "")
    looks_like_html = "html" in content_type or body.lstrip()[:14].lower().startswith(
        ("<!doctype", "<html")
    )
    if looks_like_html:
        title = re.search(r"<title>(.*?)</title>", body, re.IGNORECASE | re.DOTALL)
        return (
            'HTTP %d from %s: got an HTML page ("%s"), not JSON -- '
            "this does not look like an API endpoint; check the API URL"
            % (
                response.status_code,
                response.url,
                title.group(1).strip() if title else "no title",
            )
        )
    if len(body) > 500:
        body = body[:500] + "... (truncated)"
    return "HTTP %d from %s: %s" % (response.status_code, response.url, body)


def api_get(url, token):
    """GET from the LIMS API, exiting with a clear message on connection problems."""
    try:
        return requests.get(url, headers={"Authorization": "Token %s" % token})
    except requests.exceptions.ConnectionError as err:
        logging.error(
            "Could not connect to LIMS at %s -- is the API online and the "
            "URL correct? (%s)",
            url,
            err,
        )
        sys.exit(1)
    except requests.exceptions.Timeout as err:
        logging.error("Request to %s timed out (%s)", url, err)
        sys.exit(1)
    except requests.exceptions.RequestException as err:
        logging.error("Request to %s failed: %s", url, err)
        sys.exit(1)


def parse_json(response):
    """Parse a JSON response body, exiting with a clear message if it isn't JSON."""
    try:
        return response.json()
    except ValueError:
        logging.error(
            "Expected JSON but got something else: %s", describe_response(response)
        )
        sys.exit(1)


def get_processing_info_project(api_url, token, id, outfile):
    # get all LNs
    # then get all AGGs
    # then get all AGG info

    logging.info("Setting up project #%s", id)

    info = api_get(
        "%s/aggregation/file_detail/?library__sample__tissue_culture__project=%s&page_size=1000"
        % (api_url, id),
        token,
    )

    if info.ok:
        result = parse_json(info)
        with open(outfile, "w") as output:
            json.dump(result, output, sort_keys=True, indent=4, separators=(",", ": "))
    else:
        logging.error(
            "Could not get processing info for project %s: %s",
            id,
            describe_response(info),
        )
        sys.exit(1)

    return


def get_processing_info_experiment(api_url, token, id, outfile):
    logging.info("Setting up experiment #%s", id)

    info = api_get("%s/experiment/%s/schema" % (api_url, id), token)

    if info.ok:
        result = parse_json(info)
        with open(outfile, "w") as output:
            json.dump(result, output, indent=4, separators=(",", ": "))
    else:
        logging.error(
            "Could not get processing info for experiment %s: %s",
            id,
            describe_response(info),
        )
        sys.exit(1)

    return


def get_processing_info_alignment_group(api_url, token, group_id, outfile):
    info = api_get(
        "%s/flowcell_lane_alignment_group/%d/processing_information/"
        % (api_url, group_id),
        token,
    )

    if info.ok:
        result = parse_json(info)
        logging.info("Writing results to %s", outfile)
        with open(outfile, "w") as output:
            json.dump(result, output, sort_keys=True, indent=4, separators=(",", ": "))
    else:
        logging.error(
            "Could not get processing info for alignment group %s: %s",
            group_id,
            describe_response(info),
        )
        sys.exit(1)

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

    api_url = poptions.base_api_url or os.environ.get("LIMS_API_URL")
    if not api_url:
        logging.error(
            "Could not find LIMS API URL: pass -a/--api or set LIMS_API_URL "
            "(it may be set but empty).\n"
        )
        sys.exit(1)
    if not api_url.lower().startswith(("http://", "https://")):
        logging.error(
            'LIMS API URL "%s" does not look like a URL '
            "(expected it to start with http:// or https://)",
            api_url,
        )
        sys.exit(1)

    token = poptions.token or os.environ.get("LIMS_API_TOKEN")
    if not token:
        logging.error(
            "Could not find LIMS API TOKEN: pass -t/--token or set LIMS_API_TOKEN "
            "(it may be set but empty).\n"
        )
        sys.exit(1)

    if poptions.project:
        logging.info(
            "Getting aggregation information for project #%s", poptions.project
        )
        get_processing_info_project(api_url, token, poptions.project, poptions.outfile)

    if poptions.flowcell:
        logging.info("Getting alignment groups for %s", poptions.flowcell)

        alignment_groups = api_get(
            "%s/flowcell_lane_alignment_group/?flowcell__label=%s"
            % (api_url, poptions.flowcell),
            token,
        )

        if not alignment_groups.ok:
            logging.error(
                "Could not get alignment groups for flowcell %s: %s",
                poptions.flowcell,
                describe_response(alignment_groups),
            )
            sys.exit(1)

        results = parse_json(alignment_groups)
        if results["count"] == 0:
            # Distinguish a missing flowcell from one with no alignment groups
            flowcell_check = api_get(
                "%s/flowcell_run/?label=%s" % (api_url, poptions.flowcell),
                token,
            )
            if not flowcell_check.ok:
                logging.error(
                    "Could not find an alignment group for flowcell %s "
                    "(checking whether the flowcell exists also failed: %s)",
                    poptions.flowcell,
                    describe_response(flowcell_check),
                )
            else:
                flowcells = parse_json(flowcell_check)
                if flowcells["count"] == 0:
                    logging.error(
                        'Flowcell with label "%s" does not exist in LIMS',
                        poptions.flowcell,
                    )
                else:
                    logging.error(
                        'Flowcell "%s" exists in LIMS (id %s) '
                        "but has no alignment groups",
                        poptions.flowcell,
                        flowcells["results"][0]["id"],
                    )
            sys.exit(1)
        if results["count"] > 1:
            logging.error(
                "More than one alignment group found: %s",
                ", ".join(["%d" % ag["id"] for ag in results["results"]]),
            )
            sys.exit(1)

        get_processing_info_alignment_group(
            api_url, token, results["results"][0]["id"], poptions.outfile
        )

    if poptions.alignment_group:
        get_processing_info_alignment_group(
            api_url, token, poptions.alignment_group, poptions.outfile
        )

    if poptions.experiment:
        logging.info(
            "Getting aggregation information for experiment #%s", poptions.experiment
        )
        get_processing_info_experiment(
            api_url, token, poptions.experiment, poptions.outfile
        )


############

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
