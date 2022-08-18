#pylint disable=invalid-whitespace, invalid-name


import re
import csv
import argparse
import datetime
import hashlib
import json
import logging
import os
import sys
import time
from collections import defaultdict

sys.path.insert(
    1, os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "lims",
        "stamlims_api"
))

from stamlims_api.lims import aggregations, content_types
from stamlims_api import rest

lane_tags = None
flowcell_lane_cache = dict()
flowcell_contenttype = None

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger('upload_data.py')

script_options = {
    "base_api_url": None,
    "basedir": os.getcwd(),
    "quiet": False,
    "debug": False,

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
        help="Your authentication token.")

    parser.add_argument("sample_config",
        help="The sample_config.tsv file")
    parser.add_argument("processing_json",
        help="The processing.json file")
    parser.add_argument("--output_file_directory", default=".")


    parser.add_argument("--skip_md5", dest="skip_md5", action="store_true",
        help="Don't calculate md5sum")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


def md5sum_file(path):
    md5sum = hashlib.md5()

    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''):
            md5sum.update(chunk)

    return md5sum.hexdigest()

def url_join(*args):
    url = "/".join([ x.rstrip('/') for x in args ])
    return url

class UploadLIMS(object):

    def __init__(self, api_url, token):
        self.count_types = {}
        self.flowcelllane_contenttype = None
        self.alignment_contenttype = None
        self.aggregation_contenttype = None
        self.flowcell_lane_cache = {}
        self.api = rest.setup_api({rest.LIMS_URL_OPT_VAR: api_url,
                                   rest.LIMS_TOKEN_OPT_VAR: token,
                                   rest.RAISE_ON_ERROR_VAR: True})
        self.get_cache = {}

    def get(self, url):
        if url not in self.get_cache:
            self.get_cache[url] = self.api.get_single_result(url)
        return self.get_cache[url]

    def get_by_full_url(self, url):
        if url not in self.get_cache:
            self.get_cache[url] = self.api.get_single_result(url=url)
        return self.get_cache[url]

    def get_by_id(self, base_url, id, message=None):
        url = "%s/%d/" % (base_url, id)
        result = self.get(url)
        if not result:
            if message is None:
                message = "Failed to fetch %s" % url
            log.critical(message)
        return result

    def get_single_result(self, fetch_url, query=None, field=None):
        """
        Using a list API url that should bring up a single item, retrieve that single item if it exists.
        """
        result = self.api.get_single_list_result(url_addition=fetch_url, query_arguments=query)
        if result is None:
            return None
        if field is not None:
            return result[field]
        return result

    def get_list_result(self, url, query=None):
        return self.api.get_list_result(
            url_addition=url,
            query_arguments=query,
            item_limit=1000000,
            page_size=1000,
        )

    def put(self, *args, **kwargs):
        # TODO: s/patch/put/
        return self.api.patch_single_result(*args, **kwargs)

    def post(self, *args, **kwargs):
        return self.api.post_single_result(*args, **kwargs)

    def patch(self, *args, **kwargs):
        return self.api.patch_single_result(*args, **kwargs)

    def get_flowcell_url_by_label(self, label):
        return self.get_single_result('flowcell_run/',
                                      field = 'url',
                                      query={"label":label})

    def get_contenttype(self, contenttype_name):
        """
        Appname uses capitalization, modelname does not.
        """

        (appname, modelname) = contenttype_name.split(".")

        query = {
            'app_label': appname,
            'model': modelname,
        }
        ct = self.get_single_result('content_type/', query=query)
        if not ct:
            log.critical("Could not fetch content type %s" % contenttype_name)

        return ct

    def get_file_purpose_url(self, slug):
        return self.get_single_result('file_purpose/',
                                      query={"slug": slug},
                                      field="url")

    def get_file_type(self, slug):
        return self.get_single_result('file_type/',
                                      field="url",
                                      query={"slug":slug})


    def upload_directory_attachment(self, path, contenttype_name, object_id, file_purpose=None):
        path = os.path.abspath(path)

        if not (contenttype_name and object_id):
            log.error("Cannot attach file %s without both content type and object_id" % path)
            return False

        contenttype = self.get_contenttype(contenttype_name)

        if not contenttype:
            log.error("Cannot attach file %s without contenttype result" % path)
            return False

        purpose = self.get_file_purpose_url(file_purpose)

        if file_purpose and not purpose:
            log.error("Could not find file purpose %s for uploading directory %s" % (file_purpose, path))
            return False
        elif purpose:
            log.debug("File purpose: %s" % purpose)

        exists = self.get_single_result('directory/', query={"path":path})

        if exists:
            data = exists
        else:
            data = {}

        data.update({
            'path': path,
            'content_type': contenttype['url'],
            'object_id': object_id,
            'purpose': purpose
        })

        if exists:
            log.info("Updating information for directory %s" % path)
            result = self.put(url=data['url'], data=data)
        else:
            log.info("Uploading information for directory %s" % path)
            result = self.post("directory/", data=data)

        if not result:
            log.error("Could not upload directory %s" % path)
            log.debug(data)
        else:
            log.debug(result)

        return True

    def upload_file(self, path, contenttype_name, object_ids, file_purpose=None, file_type=None, skip_md5=False):
        log.info("Gathering data...")
        upload_data = self.get_file_upload_data(path, contenttype_name, file_purpose, file_type, skip_md5)
        if skip_md5:
            log.info("Skipping md5sum")
            upload_data['md5sum'] = '0'
        else:
            log.info("Running md5sum...")
            upload_data['md5sum'] = md5sum_file(path)

        content_type_id = re.search("(\d+)/?$", upload_data['content_type']).group(1)
        purpose_id = re.search("(\d+)/?$", upload_data['purpose']).group(1)
        for object_id in object_ids:
            data = {"object_id": object_id, **upload_data}
            exists = self.get_single_result("file/",
                    query={"object_id": object_id,
                        "purpose": purpose_id,
                        "content_type": content_type_id})

            if exists:
                log.info("Updating information for file %s: lane %d" % (path, object_id))
                result = self.put(url=exists['url'], data=data)
            else:
                log.info("Uploading information for file %s: lane %d" % (path, object_id))
                result = self.post("file/", data=data)
        
            if not result:
                log.error("Could not upload file %s for ID %d" % (path, object_id))
                log.debug(data)
            else:
                log.debug(result)



    def get_file_upload_data(self, path, contenttype_name,  file_purpose=None, file_type=None, skip_md5_check=False):
        path = os.path.abspath(path)


        contenttype = self.get_contenttype(contenttype_name)

        if not contenttype:
            log.error("Cannot attach file %s without contenttype result" % path)
            return False

        purpose = self.get_file_purpose_url(file_purpose)

        if file_purpose and not purpose:
            log.error("Could not find file purpose %s for uploading file %s" % (file_purpose, path))
            return False
        elif purpose:
            log.debug("File Purpose: %s" % purpose)

        ftype = self.get_file_type(file_type)

        if file_type and not ftype:
            log.error("Could not find file type %s for uploading file %s" % (file_type, path))
            return False
        elif purpose:
            log.debug("File Type: %s" % ftype)


        file_size = os.path.getsize(path)
        last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(path))

        #if exists:
        #recorded_mtime = datetime.datetime.fromtimestamp(time.mktime(time.strptime( exists["file_last_modified"], "%Y-%m-%dT%H:%M:%S")))

        # TODO: Make time-checking work!
        # Current issue: sub-second precision.
        data = {
            'path': path,
            'content_type': contenttype["url"],
            'purpose': purpose,
            'filetype': ftype,
            'file_last_modified': last_modified,
            'size_bytes': file_size,
        }

        log.debug(data)
        return data


    def get_flowcelllane_contenttype(self):
        if not self.flowcelllane_contenttype:
            self.flowcelllane_contenttype = self.get_contenttype('SequencingData.flowcelllane')
        return self.flowcelllane_contenttype

    def get_flowcell_lane(self, flowcell_lane_id):
        return self.get_by_id('flowcell_lane', flowcell_lane_id)

    def get_library(self, library_id):
        return self.get_by_id('library', library_id)


    def upload_altseq_flowcell(self, sample_config, processing_dict, outdir):
        # (Filepath, purpose) -> [lane_ids]
        files_to_upload = defaultdict(list)
        for row in sample_config:
            idx = row['barcode_index']
            lane = int(row['lane'])
            name = row['name']
            # Get lane IDs for each file
            lane_ids = [
                    l['id']
                    for l in processing_dict['libraries']
                    if l['barcode1']['reverse_sequence'] == idx and int(l['lane']) == lane
                    ]
            r1_file = os.path.join(outdir, "%s_R1.fq.gz" % name)
            r2_file = os.path.join(outdir, "%s_R2.fq.gz" % name)
            if not os.path.exists(r1_file):
                raise Exception("No file %s" % r1_file)
            if not os.path.exists(r2_file):
                raise Exception("No file %s" % r2_file)

            files_to_upload[(r1_file, "r1-fastq")].extend(lane_ids)
            files_to_upload[(r2_file, "r2-fastq")].extend(lane_ids)

        for ((path, purpose), lane_ids) in files_to_upload.items():
            print(path, purpose, len(lane_ids))

            self.upload_file(path,
                    "SequencingData.flowcelllane",
                    lane_ids,
                    file_purpose=purpose,
                    file_type="fastq",
                    skip_md5=True)



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

    uploader = UploadLIMS(api_url, token)

    with open(poptions.sample_config) as f:
        sample_config = [row for row in csv.DictReader(f, delimiter="\t")]
    with open(poptions.processing_json) as f:
        processing = json.loads(f.read())
    uploader.upload_altseq_flowcell(sample_config, processing, poptions.output_file_directory)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
