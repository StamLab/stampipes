#!/usr/bin/env python3
"""
Uploads all the results of alt-seq processing to LIMS
"""

import argparse
import csv
import datetime
import hashlib
import json
import logging
import os
import re
import sys
from collections import defaultdict
from functools import lru_cache

# Make sure we can load our vendored stamlims_api dependency
sys.path.insert(
    1,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "lims", "stamlims_api"
    ),
)


from stamlims_api import rest  # pylint: disable=wrong-import-position,import-error

JSON_REPORT_CLASS_SLUG = "altseq-flowcell-report-starsolo"

LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG = logging.getLogger("upload_data.py")

script_options = {
    "base_api_url": None,
    "quiet": False,
    "debug": False,
    "dry_run": False,
}


class HashableDict(dict):
    """
    A simple hashable dict
    Helps cache our GET requests even w/ query params
    """

    def __hash__(self):
        return hash(frozenset(self.items()))


def parser_setup():
    """Command-line argument setup"""
    parser = argparse.ArgumentParser()

    run_opts = parser.add_argument_group("core params")
    log_opts = parser.add_argument_group("logging options")
    lims_opts = parser.add_argument_group("lims options")

    log_opts.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Don't print info messages (only WARN and higher).",
    )
    log_opts.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="Print all debug messages.",
    )

    lims_opts.add_argument(
        "-a",
        "--api",
        dest="base_api_url",
        help="The base API url, if not the default live LIMS.",
    )
    lims_opts.add_argument(
        "-t", "--token", dest="token", help="Your authentication token."
    )

    run_opts.add_argument("sample_config", help="The sample_config.tsv file")
    run_opts.add_argument("processing_json", help="The processing.json file")
    run_opts.add_argument(
        "--output_file_directory",
        default=".",
        help="The output directory files are stored in. Defaults to cwd.",
    )

    run_opts.add_argument(
        "--skip_md5",
        dest="skip_md5",
        action="store_true",
        help="Don't calculate md5sum (debug/dev only)",
    )

    run_opts.add_argument(
        "-n",
        "--dry_run",
        dest="dry_run",
        action="store_true",
        help="Do not upload anything to LIMS, instead print actions that would be taken",
    )

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


def md5sum_file(path):
    """Calculates the md5sum of a file's contents"""
    md5sum = hashlib.md5()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5sum.update(chunk)

    return md5sum.hexdigest()


def parse_counts_file(counts_file: str):
    """
    Given a file name, reads a stats file
    format: one stat per line: `name value` (separated by whitespace)
    returns a dict of str->int
    """
    stats = {}
    with open(counts_file, "r") as counts:
        for line in counts:
            values = line.split()
            count_type_name = values[0]
            if not count_type_name:
                continue
            count = int(values[1])
            stats[count_type_name] = count
    return stats


def build_counts(alignment_id, counts_file):
    """
    Convert stats into a form ready to be uploaded to LIMS with the
    bulk-stat-create endpoint
    """
    parsed_stats = parse_counts_file(counts_file)
    return {
        "object_id": alignment_id,
        "content_type": "SequencingData.flowcelllanealignment",
        "stats": parsed_stats,
    }


class UploadLIMS:
    """
    Contains the logic for uploading things to LIMS
    Uses caching for most GET requests
    """

    def __init__(self, api_url, token, dry_run=False, skip_md5=False):
        # self.count_types = {}
        # self.flowcelllane_contenttype = None
        # self.alignment_contenttype = None
        # self.aggregation_contenttype = None
        self.api = rest.setup_api(
            {
                rest.LIMS_URL_OPT_VAR: api_url,
                rest.LIMS_TOKEN_OPT_VAR: token,
                rest.RAISE_ON_ERROR_VAR: True,
            }
        )
        self.dry_run = dry_run
        self.skip_md5 = skip_md5

    @lru_cache(maxsize=None)
    def get(self, url):
        """Cached version of api.get_single_result"""
        return self.api.get_single_result(url)

    def get_by_id(self, base_url, object_id, err_message=None):
        """Constructs url from ID and calls get"""
        url = "%s/%d/" % (base_url, object_id)
        result = self.get(url)
        if not result:
            if err_message is None:
                err_message = "Failed to fetch %s" % url
            LOG.critical(err_message)
        return result

    @lru_cache(maxsize=None)
    def _get_single_result(self, fetch_url, query=None, field=None):
        """Internal memo-izable function, do not use directly"""
        result = self.api.get_single_list_result(
            url_addition=fetch_url, query_arguments=query
        )
        if result is None:
            return None
        if field is not None:
            return result[field]
        return result

    def get_single_result(self, fetch_url, query=None, field=None):
        """
        Using a list API url that should bring up a single item, retrieve that
        single item if it exists.
        """
        if isinstance(query, dict) and not isinstance(query, HashableDict):
            query = HashableDict(query)
        return self._get_single_result(fetch_url, query, field)

    # Not currently used
    @lru_cache(maxsize=None)
    def _get_list_result(self, url, query=None):
        return self.api.get_list_result(
            url_addition=url,
            query_arguments=query,
            item_limit=1000000,
            page_size=1000,
        )

    def get_list_result(self, url, query=None):
        if isinstance(query, dict) and not isinstance(query, HashableDict):
            query = HashableDict(query)
            LOG.debug("Query is now: %s", query)
        return self._get_list_result(url, query)

    def put(self, *args, **kwargs):
        """
        PUT data to LIMS
        """
        if self.dry_run:
            LOG.info("Dry run, would have put %s, %s", args, kwargs)
            return None
        # FIXME: Should use PUT method once API lib supports it
        return self.api.patch_single_result(*args, **kwargs)

    def post(self, *args, **kwargs):
        """
        POST data to LIMS
        """
        if self.dry_run:
            LOG.info("Dry run, would have post %s, %s", args, kwargs)
            return None
        return self.api.post_single_result(*args, **kwargs)

    def patch(self, *args, **kwargs):
        if self.dry_run:
            LOG.info("Dry run, would have patch %s, %s", args, kwargs)
            return None
        return self.api.patch_single_result(*args, **kwargs)

    # def get_flowcell_url_by_label(self, label):
    #     return self.get_single_result(
    #         "flowcell_run/", field="url", query={"label": label}
    #     )

    def get_contenttype(self, contenttype_name):
        """
        Appname uses capitalization, modelname does not.
        """

        (appname, modelname) = contenttype_name.split(".")

        query = {
            "app_label": appname,
            "model": modelname,
        }
        ct = self.get_single_result("content_type/", query=query)
        if not ct:
            LOG.critical("Could not fetch content type %s", contenttype_name)

        return ct

    def get_file_purpose_url(self, slug):
        """Get file purpose url from slug"""
        return self.get_single_result(
            "file_purpose/", query={"slug": slug}, field="url"
        )

    def get_file_type_url(self, slug):
        """Gets the file type URL for a slug"""
        return self.get_single_result("file_type/", field="url", query={"slug": slug})

    def upload_directory_attachment(
        self, path, contenttype_name, object_id, file_purpose=None
    ):
        """Uploads a single directory to a LIMS object"""
        path = os.path.abspath(path)
        if not (contenttype_name and object_id):
            LOG.error(
                "Cannot attach file %s without both content type and object_id", path
            )
            return False

        contenttype = self.get_contenttype(contenttype_name)
        if not contenttype:
            LOG.error("Cannot attach file %s without contenttype result", path)
            return False

        purpose = self.get_file_purpose_url(file_purpose)
        if file_purpose and not purpose:
            LOG.error(
                "Could not find file purpose %s for uploading directory %s",
                file_purpose,
                path,
            )
            return False
        LOG.debug("File purpose: %s", purpose)

        existing_data = self.get_single_result("directory/", query={"path": path})
        data = existing_data if existing_data else {}

        data.update(
            {
                "path": path,
                "content_type": contenttype["url"],
                "object_id": object_id,
                "purpose": purpose,
            }
        )

        if existing_data:
            LOG.info("Updating information for directory %s", path)
            result = self.put(url=data["url"], data=data)
        else:
            LOG.info("Uploading information for directory %s", path)
            result = self.post("directory/", data=data)

        if not result:
            LOG.error("Could not upload directory %s", path)
            LOG.debug(data)
        else:
            LOG.debug(result)

        return True

    def upload_file(
        self, path, contenttype_name, object_ids, file_purpose=None, file_type=None
    ):
        """
        Upload a file's metadata to LIMS
        It will be attached to many objects.
        """
        # FIXME: This method makes a GET and PUT request for every single object
        # Will require LIMS API updates to enable a more performant solution

        upload_data = self.get_file_upload_data(
            path, contenttype_name, file_purpose, file_type
        )
        LOG.debug("Uploading file %s, to %d objects", path, len(object_ids))
        if self.skip_md5:
            LOG.info("Skipping md5sum")
            upload_data["md5sum"] = "0"
        else:
            LOG.debug("Running md5sum...")
            upload_data["md5sum"] = md5sum_file(path)

        content_type_id = re.search(r"(\d+)/?$", upload_data["content_type"]).group(1)
        purpose_id = re.search(r"(\d+)/?$", upload_data["purpose"]).group(1)
        for object_id in object_ids:
            upload_data.update({"object_id": object_id})
            exists = self.get_single_result(
                "file/",
                query={
                    "object_id": object_id,
                    "purpose": purpose_id,
                    "content_type": content_type_id,
                },
            )

            if exists:
                if exists == upload_data:
                    LOG.info(
                        "No change to information for file %s, lane %d, not updating",
                        path,
                        object_id,
                    )
                    result = True
                else:
                    LOG.info(
                        "Updating information for file %s: lane %d", path, object_id
                    )
                    result = self.put(url=exists["url"], data=upload_data)
            else:
                LOG.info(
                    "Uploading information for file %s: lane %d, data=%s",
                    path,
                    object_id,
                    upload_data,
                )
                result = self.post("file/", data=upload_data)

            if not result:
                LOG.error("Could not upload file %s for ID %d", path, object_id)
                LOG.debug(upload_data)
            else:
                LOG.debug(result)

    def get_file_upload_data(
        self, path, contenttype_name, file_purpose=None, file_type=None
    ):
        """
        Gets the file upload data that is easy to query
        (notable omission: md5sum, as it takes a long time to calculate)
        """
        path = os.path.abspath(path)

        contenttype = self.get_contenttype(contenttype_name)
        if not contenttype:
            LOG.error("Cannot attach file %s without contenttype result", path)
            return False

        purpose = self.get_file_purpose_url(file_purpose)
        if file_purpose and not purpose:
            LOG.error(
                "Could not find file purpose %s for uploading file %s",
                file_purpose,
                path,
            )
            return False
        if purpose:
            LOG.debug("File Purpose: %s", purpose)

        ftype = self.get_file_type_url(file_type)
        if file_type and not ftype:
            LOG.error(
                "Could not find file type %s for uploading file %s", file_type, path
            )
            return False
        if file_type:
            LOG.debug("File Type: %s", ftype)

        file_size = os.path.getsize(path)
        last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(path))

        # Current issue: sub-second precision.
        data = {
            "path": path,
            "content_type": contenttype["url"],
            "purpose": purpose,
            "filetype": ftype,
            "file_last_modified": last_modified,
            "size_bytes": file_size,
        }

        LOG.debug(data)
        return data

    def get_flowcell_lane(self, flowcell_lane_id):
        """Gets the flowcell lane by ID"""
        return self.get_by_id("flowcell_lane", flowcell_lane_id)

    def get_library(self, library_id):
        """Gets the library by ID (NOT library number)"""
        return self.get_by_id("library", library_id)

    def upload_flowcell_report(self, data):
        flowcell_labels = set(pool["flowcell_label"] for pool in data)
        assert len(flowcell_labels) == 1
        flowcell_label = flowcell_labels.pop()

        report_name = "Alt-seq stats: FC%s" % flowcell_label

        flowcell_lims_info = self.get_single_result(
            "flowcell_run/?label=%s" % flowcell_label
        )
        content_type_id = flowcell_lims_info["object_content_type"]
        content_type = self.get_by_id("content_type", content_type_id)
        object_id = flowcell_lims_info["id"]
        json_report_class = self.get_single_result(
            "json_report_class/", query={"slug": JSON_REPORT_CLASS_SLUG}
        )

        # See if report already exists
        existing_reports = self.get_list_result(
            "json_report/",
            query={
                "object_id": object_id,
                "content_type": content_type["id"],
                "report_class": json_report_class["id"],
                "page_size": 2,
            },
        )

        data_to_send = {
            "object_id": object_id,
            "content_type": content_type["url"],
            "report_class": json_report_class["url"],
            "name": report_name,
            "json_content": json.dumps(data),
        }
        if len(existing_reports) == 0:
            self.post("json_report/", data=data_to_send)
            # No report exists yet, upload a new one
        elif len(existing_reports) == 1:
            # Exactly one report, update it
            url_to_patch = "json_report/%d/" % existing_reports[0]["id"]
            self.patch(url_to_patch, data=data_to_send)
        else:
            # Error! too many reports
            LOG.critical("Too many JSON reports exist")
            raise "Too many JSON reports exist, exiting"

    def upload_altseq_flowcell(self, sample_config, processing_dict, outdir):
        """
        Main function for this script.
        Given paths to the sample_config file, processing_dict, and outdir,
        upload to LIMS:
        1) Paths for fastq files for each lane
        # 2) Stats for each alignment
        3) Flowcell-level pool stats
        """
        # (Filepath, purpose) -> [lane_ids]
        files_to_upload = defaultdict(list)

        # Augment processing_dict with sample_config info
        processing_info = []
        for row in sample_config:
            barcode_index = row["barcode_index"]
            lane = int(row["lane"])
            pool_name = row["pool_name"]
            sample_name = row["sample_name"]
            for idx, lib in enumerate(processing_dict["libraries"]):
                if int(lib["lane"]) == lane and lib["barcode_index"] == barcode_index:
                    lib.update({"pool_name": pool_name, "sample_name": sample_name})
                    processing_info.append(lib)

        # TODO: Doesn't yet make use of the above augmented info
        for row in sample_config:
            (idx, _otheridx) = row["barcode_index"].split("-")
            lane = int(row["lane"])
            name = row["pool_name"]
            LOG.debug("idx=%s, lane=%d, name=%s", idx, lane, name)
            # Get lane IDs for each file
            lane_ids = [
                lib["id"]
                for lib in processing_dict["libraries"]
                if lib["barcode1"]["reverse_sequence"] == idx
                and int(lib["lane"]) == lane
            ]
            r1_file = os.path.join(outdir, name, "R1.fq.gz")
            r2_file = os.path.join(outdir, name, "R2.fq.gz")
            if not os.path.exists(r1_file):
                raise Exception("No file %s" % r1_file)
            if not os.path.exists(r2_file):
                raise Exception("No file %s" % r2_file)

            files_to_upload[(r1_file, "r1-fastq")].extend(lane_ids)
            files_to_upload[(r2_file, "r2-fastq")].extend(lane_ids)

        # Upload files.
        for (path, purpose), lane_ids in files_to_upload.items():
            # print(path, purpose, len(lane_ids))
            self.upload_file(
                path,
                "SequencingData.flowcelllane",
                list(set(lane_ids)),
                file_purpose=purpose,
                file_type="fastq",
            )

        # Commented out because we aren't making alignments for these...
        # # Now upload counts.
        # # We can do this all as one call.
        # # (Assuming LIMS doesn't time out)
        # all_counts = []
        # for lib in processing_info:
        #     if not len(lib["alignments"]) == 1:
        #         LOG.critical("Lib must have exactly 1 aligment %s", lib)
        #     align_id = lib["alignments"][0]["id"]
        #     counts_file = os.path.join(
        #         outdir,
        #         lib["pool_name"],
        #         "analysis",
        #         "Gene",
        #         "%s.stats.txt" % lib["sample_name"],
        #     )
        #     all_counts.append(build_counts(align_id, counts_file))
        # # print(json.dumps(all_counts))
        # self.post("stats/create/", all_counts)

        with open(os.path.join(outdir, "flowcell_stats.json")) as json_file:
            flowcell_data = json.loads(json_file.read())
            self.upload_flowcell_report(flowcell_data)


def main():
    """
    This is the main body of the program that uses the arguments from the
    command line.
    """

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=LOG_FORMAT)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
    else:
        # Set up the default logging levels
        logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
        # Make this a little less noisy by default
        requests_log = logging.getLogger("requests.packages.urllib3.connectionpool")
        requests_log.setLevel(logging.WARN)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
        LOG.debug("Using LIMS API endpoint: %s from environment", api_url)
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
        LOG.debug("Using LIMS API endpoint: %s from options", api_url)
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

    uploader = UploadLIMS(
        api_url, token, dry_run=poptions.dry_run, skip_md5=poptions.skip_md5
    )

    with open(poptions.sample_config) as f:
        sample_config = list(csv.DictReader(f, delimiter="\t"))
    with open(poptions.processing_json) as f:
        processing = json.loads(f.read())
    uploader.upload_altseq_flowcell(
        sample_config, processing, poptions.output_file_directory
    )


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell
# after importing without automatically running it
if __name__ == "__main__":
    main()
