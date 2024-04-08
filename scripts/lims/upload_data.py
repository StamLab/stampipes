#pylint disable=invalid-whitespace, invalid-name

import argparse
import datetime
import hashlib
import json
import logging
import os
import re
import sys
import time
from zipfile import ZipFile

sys.path.insert(
    1, os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
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

    "aggregation_id": None,
    "start_aggregation": False,
    "complete_aggregation": False,
    "clear_aggregation_stats": False,

    "alignment_id": None,
    "flowcell": None,
    "flowcell_lane_id": None,

    "fastqc_counts": False,
    "fastqc_files": [],

    "spot_file": None,
    "spot_dup_file": None,
    "dups_file": None,
    "counts_file": None,
    "rna_file": None,
    "barcode_report_file": None,

    "version_file": None,
    "adapter_file": None,

    "align_start_time": False,
    "align_complete_time": False,

    "attach_file": None,
    "attach_directory": None,
    "attach_file_contenttype": None,
    "attach_file_objectid": None,
    "attach_file_purpose": None,
    "attach_file_type": None,

    "clear_align_stats": False,

    "skip_md5_check": False,
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
        help="The flowcell we're working on.  Enter it to clear cache after uploading.")

    parser.add_argument("--aggregation_id", dest="aggregation_id", type=int)
    parser.add_argument("--clear_aggregation_stats", dest="clear_aggregation_stats", action="store_true",
        help="Clear the statistics/files for a given aggregation.")
    parser.add_argument("--start_aggregation", dest="start_aggregation", action="store_true",
        help="Set the current time for the alignment's start time.")
    parser.add_argument("--complete_aggregation", dest="complete_aggregation", action="store_true",
        help="Set the current time for the alignment's complete time.")

    parser.add_argument("--clear_align_stats", dest="clear_align_stats", action="store_true",
        help="Clear the statistics/files for a given alignment.")

    # these should go together
    parser.add_argument("--alignment_id", dest="alignment_id", type=int)
    parser.add_argument("--spotfile", dest="spot_file",
        help="The SPOT output file.  Best paired with a dupfile.  Needs alignment id.")
    parser.add_argument("--spotdupfile", dest="spot_dup_file",
        help="The Picard dup results file paired with a spotfile.  Needs alignment id.")

    # requires alignment_id
    parser.add_argument("--start_alignment_progress", dest="align_start_time", action="store_true",
        help="Set the current time for the alignment's start time.")
    parser.add_argument("--finish_alignment", dest="align_complete_time", action="store_true",
        help="Set the current time for the alignment's complete time.")

    # also needs alignment_id
    parser.add_argument("--countsfile", dest="counts_file",
        help="A tab delineated list of counts.  Needs alignnment id.")

    # requires alignment_id
    parser.add_argument("--version_file", dest="version_file",
        help="A version file for alignments.")
    parser.add_argument("--adapter_file", dest="adapter_file",
        help="An adapter file for alignments.")

        # A lane can have multiple fastQC files, one for each read
    parser.add_argument("--flowcell_lane_id", dest="flowcell_lane_id", type=int,
        help="The ID of the flowcell lane we're working on.")
    parser.add_argument("--fastqcfile", dest="fastqc_files", action="append",
        help="A FastQC ZIP file to upload.")
    parser.add_argument("--insertsfile", dest="inserts_file",
        help="A Picard CollectInsertSizeMetrics text file for an alignment.")
    parser.add_argument("--dupsfile", dest="dups_file",
        help="A Picard MarkDuplicates text file for an alignment.")
    parser.add_argument("--rnafile", dest="rna_file",
        help="The RNA metric output file")
    parser.add_argument("--fastqc_counts", dest="fastqc_counts", action="store_true",
        help="Use the given fastqc files to create total/pf/qc counts. Must have an alignment id.")
    parser.add_argument("--barcode_report", dest="barcode_report_file",
        help="The barcode report JSON file")

    parser.add_argument("--attach_file", dest="attach_file",
        help="The full path to a file to attach to a LIMS object.")
    parser.add_argument("--attach_directory", dest="attach_directory",
        help="The full path to a directory to attach to a LIMS object.")
    parser.add_argument("--attach_file_contenttype", dest="attach_file_contenttype",
        help="The content type to attach to, aka SequencingData.flowcelllanealignment")
    parser.add_argument("--attach_file_objectid", dest="attach_file_objectid", type=int,
        help="The object ID to attach to.")
    parser.add_argument("--attach_file_purpose", dest="attach_file_purpose",
        help="The file's purpose slug.")
    parser.add_argument("--attach_file_type", dest="attach_file_type",
        help="The file's type slug.")

    parser.add_argument("--skip_md5_check", dest="skip_md5_check", action="store_true",
        help="If file exists and path/size match, don't check md5sum.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

def split_sample_name(samplename):

    m = re.match(r'(?P<sample>[^/]+)_(?P<barcode>[AGTC-]+|NoIndex)_L00(?P<lane>[0-9])', samplename)

    if not m:
        log.error("Could not parse sample name: %s" % samplename)
        return None

    return { "sample": m.group('sample'), "barcode": m.group('barcode'), "lane": m.group('lane') }

def get_spot_score(spot_file):

    contents = open(spot_file, 'r').read()
    stats = contents.split("\n")[1].split()

    return {"total_tags": int(stats[0]), "tags_in_hotspots": int(stats[1]), "spot_score": stats[2]}

def get_dup_score(spotdup_file):
    if not spotdup_file:
        return None

    infile = open(spotdup_file, 'r')

    try:
        for line in infile:
            if line.startswith("LIBRARY"):
                percent_duplication = float(next(infile).strip().split("\t")[8])

        return percent_duplication
    except UnboundLocalError as e:
        log.error("Unbound Local Error for %s" % spotdup_file)
        log.error(e)
    except IndexError as e:
        log.error(e)

    return None

def get_fastqc_counts(fastqc_input):

    total_m = re.search(r'Total Sequences\t(?P<total>\d+)', fastqc_input)

    if not total_m:
        log.error("Could not get total sequences from fastqc_input")
        return None

    filtered_m = re.search(r'Filtered Sequences\t(?P<filtered>\d+)', fastqc_input)

    if not filtered_m:
        log.error("Could not get filtered sequences from fastqc_input")
        return None

    return {
        'total': int(total_m.group('total')),
        'filtered': int(filtered_m.group('filtered')),
    }

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
        self.fastqc_tags = None
        self.count_types = {}
        self.flowcelllane_contenttype = None
        self.alignment_contenttype = None
        self.aggregation_contenttype = None
        self.flowcell_lane_cache = {}
        self.alignment_counts = {}
        self.picard_metrics = None
        self.fastqc_counts = {}
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

    def clear_flowcell_cache(self, flowcell):
        url = self.get_flowcell_url_by_label(flowcell)
        if url is None:
            log.error("Failure to reset flowcell cache for %s" % flowcell)
            return
        log.debug(self.api.post_single_result(url=url + "clear_cache/"))

    def clear_alignment_stats(self, alignment_id):
        url = "flowcell_lane_alignment/%d/clear_stats/" % alignment_id
        log.debug("Clearing stats: %s" % url)
        results = self.post(url)
        if results is None:
            log.error("Could not clear alignment stats for ALN%s" % alignment_id)

    def clear_aggregation_stats(self, aggregation_id):
        url = "aggregation/%d/clear_stats/" % aggregation_id
        log.debug("Clearing stats: %s" % url)
        results = self.post(url)
        if results is None:
           log.error("Could not clear aggregation stats for AGG%s" % aggregation_id)

    def start_aggregation(self, aggregation_id):
        url = "aggregation/%d/" % aggregation_id
        data = {
            "processing_started": datetime.datetime.now(),
        }
        results = self.patch(url, data=data)
        if results is None:
            log.error("Could not 'start' AGG%s" % aggregation_id)

    def complete_aggregation(self, aggregation_id):
        url = "aggregation/%d/" % aggregation_id

        data = {
            "processing_completed": datetime.datetime.now(),
            "needs_reprocessing": False,
        }

        results = self.patch(url, data=data)

        if results is None:
            log.error("Could not complete AGG%s" % aggregation_id)

    def get_fastqc_tags(self):

        if not self.fastqc_tags:
            tags = self.get_list_result('fastqc_tag/')
            if tags is None:
                log.critical("Could not fetch fastqc tags from LIMS")

            self.fastqc_tags = dict([(tag['slug'], tag) for tag in tags])

        return self.fastqc_tags

    def get_picard_metrics(self):

        if not self.picard_metrics:
            picard_metrics = self.get_list_result('picard_metric/')
            if picard_metrics is None:
                log.critical("Could not fetch picard metrics from LIMS")

            self.picard_metrics = dict([(metric['name'], metric) for metric in picard_metrics])

        return self.picard_metrics

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

    def upload_file_attachment(self, path, contenttype_name, object_id, file_purpose=None, file_type=None, skip_md5_check=False):
        path = os.path.abspath(path)

        log.info("Attaching file %s to object %d (contenttype %s)" % (path, object_id, contenttype_name))

        if not (contenttype_name and object_id):
            log.error("Cannot attach file %s without both content type and object_id" % path)
            return False

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

        exists = self.get_single_result("file/",
                                        query={"object_id": object_id,
                                               "purpose__slug": file_purpose,
                                               "content_type": contenttype['id']})

        file_size = os.path.getsize(path)
        last_modified = datetime.datetime.fromtimestamp(os.path.getmtime(path))

        if skip_md5_check and exists and exists["size_bytes"] == file_size:
            recorded_mtime = datetime.datetime.fromtimestamp(time.mktime(time.strptime(
                exists["file_last_modified"], "%Y-%m-%dT%H:%M:%S")))
            # Allow for sloppiness in NFS timestamps
            difference = recorded_mtime - last_modified
            if timedelta(minutes=-1) <= difference <= timedelta(minutes=1):
                log.info("File exists and matches recorded size, skipping %s" % path)
            return

        md5sum = md5sum_file(path)

        log.info("MD5sum: %s\tFile size: %d\tLast modified: %s" % (md5sum, file_size, str(last_modified)))

        data = {
            'path': path,
            'content_type': contenttype["url"],
            'object_id': object_id,
            'purpose': purpose,
            'filetype': ftype,
            'md5sum': md5sum,
            'file_last_modified': last_modified,
            'size_bytes': file_size,
        }

        log.debug(data)

        if exists:
            log.info("Updating information for file %s" % path)
            result = self.put(url=exists['url'], data=data)
        else:
            log.info("Uploading information for file %s" % path)
            result = self.post("file/", data=data)

        if not result:
            log.error("Could not upload file %s" % path)
            log.debug(data)
        else:
            log.debug(result)

    def get_flowcelllane_contenttype(self):
        if not self.flowcelllane_contenttype:
            self.flowcelllane_contenttype = self.get_contenttype('SequencingData.flowcelllane')
        return self.flowcelllane_contenttype

    def get_alignment_contenttype(self):
        self.alignment_contenttype = self.get_contenttype('SequencingData.flowcelllanealignment')
        return self.alignment_contenttype

    def get_aggregation_contenttype(self):
        self.aggregation_contenttype = self.get_contenttype('AggregationData.aggregation')
        return self.aggregation_contenttype

    def create_count_type(self, name):

        log.info("Creating count type %s" % name)

        is_mapq = name.startswith("mapq")
        is_samflag = name.startswith("samflag")
        is_alignment = "readlength" in name

        is_chromosome = not(is_mapq or is_samflag or is_alignment)

        data = {
            "is_chromosome": is_chromosome,
            "is_samflag": is_samflag,
            "is_mapq": is_mapq,
            "is_alignment": is_alignment,
            "codename": name,
            "title": name,
        }

        result = self.post("flowcell_lane_count_type/", data=data)

        if result.ok:
            self.count_types[name] = result
        else:
            self.count_types[name] = None
            log.warn("Could not create count type %s (%s)" % (name, str(result)))
        return self.count_types[name]


    # TODO : make sure that no more of one count type exists
    def get_alignment_counts(self, alignment_id):

        log.info("Getting alignment counts for %d" % alignment_id)
        if not alignment_id in self.alignment_counts:
            counts = self.get_list_result('flowcell_lane_count/',
                                          query={"alignment":alignment_id})
            if counts is None:
                log.critical("Could not get counts for ALN%d" % alignment_id)
            self.alignment_counts[alignment_id] = dict([(count['count_type_name'], count) for count in counts])

        return self.alignment_counts[alignment_id]

    def get_flowcell_lane(self, flowcell_lane_id):
        return self.get_by_id('flowcell_lane', flowcell_lane_id)

    def get_library(self, library_id):
        return self.get_by_id('library', library_id)

    def get_aggregation(self, aggregation_id):
        return self.get_by_id('aggregation', aggregation_id)

    def get_rna_metrics(self, alignment_id):
        exists = self.get_single_result('rna_alignment_metrics/', query={"alignment": alignment_id})
        if not exists:
            log.error("Error finding RNA metrics for alignment %d" % alignment_id)
        return exists

    def upload_rna_metrics(self, alignment_id, rna_file):
        content = open(rna_file, 'r')
        metrics = dict()
        for line in content:
            values = line.split()
            metric_name = values[0]
            metric_value = values[1]
            metrics[metric_name] = metric_value
        content.close()

        exists = self.get_rna_metrics(alignment_id)

        if exists:
            data = exists
        else:
            data = {}

        data.update({
            "alignment":               "%s/flowcell_lane_alignment/%d/" % (self.api.api_url, alignment_id),
            "input_reads":             metrics[r'input_reads'],
            "mapped_reads":            metrics[r'mapped'],
            "percent_rRNA":            metrics[r'%rRNA'],
            "percent_duplicates":      metrics[r'%duplicates'],
            "exon_intron":             metrics[r'exon:intron'],
            "percent_intergenic":      metrics[r'%intergenic'],
            "percent_chrM":            metrics[r'%chrM'],
            "percent_correct_strand":  metrics[r'%correct_strand']
        })

        if exists:
            # Currently (2014-12-22) this will fail, but that's a TODO on the LIMS side.
            log.info("Updating RNA metrics for alignment ID %d" % alignment_id)
            result = self.put(url=data['url'], data=data)
        else:
            log.info("Uploading RNA metrics for alignment ID %d" % alignment_id)
            result = self.post("rna_alignment_metrics/", data=data)
        log.debug(result)
        if not result:
            log.error("Could not upload RNA stats")

    def upload_barcode_report(self, barcode_file):
        datastring = open(barcode_file, 'r').read()
        try:
            jsondata = json.loads(datastring)
        except ValueError:
            log.error("Barcode report %s is not valid JSON" % barcode_file)
            return

        if jsondata['Sequencer'] == 'MiniSeq':
            print(jsondata['BaseDir'])
            flowcell_label = re.search( '.*_[AB](000[A-Z0-9]{6}).*$', jsondata['BaseDir'] ).group(1)
            print(flowcell_label)

        else:
            # make this more flexible eventually
            flowcell_label = re.search( '.*_[AB]([A-Z0-9]{9})$', jsondata['BaseDir'] ).group(1)

        flowcell_url = self.get_flowcell_url_by_label(flowcell_label)

        data = {
            "flowcell": flowcell_url,
            "json_data": datastring
        }

        # TODO: Don't upload redundant barcodes.
        result = self.post("barcode_report/", data=data)
        log.debug(result)

    def bulk_upload_counts(self, alignment_id, stats):
        # TODO: This isn't ready yet.
        data = {
            "object_id": alignment_id,
            "content_type": "flowcelllanealignment",
            "stats": stats,
        }
        response = self.api.post_single_result(url_addition="stat/create", json=data)
        return response

    def upload_counts(self, alignment_id, counts_file):
        parsed = self.parse_counts(counts_file)
        response = self.bulk_upload_counts(alignment_id, self.parse_counts(counts_file))
        if response is None:
            log.error("Bulk upload failed: Counts file {} for ALN{}".format(counts_file, alignment_id))
        else:
            log.info("Upload successful.")
        return
        # TODO: Remove below code

        #log.warn("Counts: %s", self.get_list_result(
        #        'flowcell_lane_count/', query={"alignment":alignment_id}
        #))
        existing_counts = {
            count['count_type_name']: (count['count'], count['url'])
            for count in self.get_list_result(
                'flowcell_lane_count/', query={"alignment":alignment_id}
            )
        }
        #log.warn("Count types: %s", self.get_list_result("flowcell_lane_count_type"))

        lane_count_types = {
            ct['codename']: ct['url']
            for ct in self.get_list_result("flowcell_lane_count_type")
        }

        for (key, value) in parsed.items():
            try:
                data = {
                    "alignment": "%s/flowcell_lane_alignment/%d/" % (self.api.api_url, alignment_id),
                    "count_type": lane_count_types[key],
                    "count": value,
                }
            except KeyError:
                log.error("Count type %s not in LIMS", key)

            existing_value = existing_counts.get(key, None)

            if existing_value is None:
                self.post(
                    "flowcell_lane_count/",
                    json=data,
                )
            elif value != existing_value[0]:
                self.patch(
                    url=existing_value[1],
                    json=data,
                )
            # else we don't need to do anything


    def parse_counts(self, counts_file):
        stats = {}
        with open(counts_file, 'r') as counts:
            for line in counts:
                values = line.split()
                count_type_name = values[0]
                if not count_type_name:
                    continue
                count = int(values[1])
                stats[count_type_name] = count
        return stats

    def upload_alignment_records(self, alignment_id, adapter_file=None, version_file=None, start_time = False, complete_time = False):

        log.info("Uploading alignment records for %d" % alignment_id)

        if not (adapter_file or version_file or start_time or complete_time):
            log.debug("No data to upload.")
            return False

        alignment = self.get_by_id("flowcell_lane_alignment", alignment_id)

        if version_file:
            alignment["versions"] = open(version_file, 'r').read()

        if adapter_file:
            alignment["trim_adapters"] = open(adapter_file, 'r').read()

        if start_time:
            alignment["start_time"] = datetime.datetime.now()

        if complete_time:
            alignment["complete_time"] = datetime.datetime.now()

        result = self.patch(url=alignment['url'], data=alignment)

        if result:
            log.info("Alignment %d updated" % alignment_id)
            log.debug(result)
        else:
            log.debug("No result for uploading %s to %s" % (str(alignment), alignment['url']))

        return True

    def upload_spot(self, alignment_id, spot_file, dup_file):

        if not spot_file and dup_file:
            log.error("Error, do not have both files for alignment %s" % alignment_id)

        spot_stats = get_spot_score(spot_file)
        percent_dup = get_dup_score(dup_file)

        data = {
            "alignment": "%s/flowcell_lane_alignment/%d/" % (self.api.api_url, alignment_id)
        }

        if spot_stats:
            data.update(spot_stats)

        data["percent_duplication"] = percent_dup

        log.debug(data["percent_duplication"])

        origspots = self.get_list_result("flowcell_lane_spot/",
                                          query={"alignment": alignment_id})
        if len(origspots) > 1:
            log.error("Could not figure out which SPOT score to upload to!")
        elif len(origspots) == 0:
            log.info("Uploading new spot for %d" % alignment_id)
            result = self.post("flowcell_lane_spot/", data=data)
            if not result:
                log.error("Could not upload SPOT")
        else:
            origspot = origspots[0]
            if (data["spot_score"] != origspot["spot_score"]
                or data["total_tags"] != origspot["total_tags"]
                or data["tags_in_hotspots"] != origspot["tags_in_hotspots"]
                or data["percent_duplication"] != origspot["percent_duplication"]
            ):
                log.info("Updating SPOT score for %d" % alignment_id)
                result = self.patch(url=origspot['url'], data=data)
                if not result:
                    log.error("Could not upload SPOT")

    def get_fastqc_contents(self, filename):

        file_in_zip = "%s/fastqc_data.txt" % os.path.splitext(os.path.basename(filename))[0]
        with ZipFile(filename) as fastqc_zip:
            with fastqc_zip.open(file_in_zip) as fastqc_report:
               return fastqc_report.read()

        return None

    def upload_fastqc(self, flowcell_lane_id, filename):

        if not self.fastqc_tags:
            self.fastqc_tags = self.get_fastqc_tags()
        if not self.flowcelllane_contenttype:
            self.flowcelllane_contenttype = self.get_flowcelllane_contenttype()

        m = re.search(r'(?P<samplename>[^/]+)_(?P<barcode>[AGTC-]+|NoIndex)_L00(?P<lane>[0-9])_(?P<read>R[12])', filename)

        if not m:
            log.error("Could not figure out information for %s" % filename)
            return False

        log.info(m.groups())

        fastqc_report = self.get_fastqc_contents(filename)

        if not fastqc_report:
           log.error("Could not read fastqc report %s" % filename)
           return False

        samplename = m.group('samplename')
        read = m.group('read')

        lane_info = self.get_flowcell_lane(flowcell_lane_id)

        if not lane_info:
            return False

        tag = self.fastqc_tags[read]

        upload = dict()

        upload['tags'] = [tag['url']]
        upload['raw_data'] = fastqc_report
        upload['content_type'] = self.flowcelllane_contenttype["url"]
        upload['object_id'] = lane_info['id']
        upload['label'] = "FC%s %s %s %s %s" % (lane_info['flowcell_label'], samplename, str(lane_info["lane"]), lane_info["barcode_index"], read)

        # does this report already exist?
        report = self.get_single_result(
            'fastqc_report/',
            query={
                "label": upload['label'],
                "object_id": upload['object_id'],
                "content_type": self.get_flowcelllane_contenttype()['id']
            })

        if report:
            # replace content
            if 'raw_data' not in report or report['raw_data'] != upload['raw_data']:
                log.info("Updating report %s" % upload['label'])
                result = self.patch(url=report['url'], data=upload)

                if result:
                    log.debug(result)
                else:
                    log.error("Could not update FastQC report %s" % report['url'])
        else:
            log.info("Uploading new fastqc report %s" % upload['label'])
            result = self.post("fastqc_report/", data=upload)

            if result:
                log.debug(result)
            else:
                log.error("Could not upload new FastQC report")

    def upload_fastqc_counts(self, alignment_id):

        if not alignment_id:
            logging.critical("Could not upload fastqc_counts without an alignment id given")
            return

        self.get_alignment_counts(alignment_id)

        total = 0
        filtered = 0

        for fastqc_file, fastqc_counts in self.fastqc_counts.items():

            if not fastqc_counts:
                log.error("Could not get counts from %s for uploading" % fastqc_file)
                return

            total += fastqc_counts["total"]
            filtered += fastqc_counts["filtered"]

        # FastQC's definition of total differs from ours
        counts = {
            "total": total + filtered,
            "qc": filtered,
            "pf": total
        }

        if not self.bulk_upload_counts(alignment_id, counts):
            log.error("Could not upload FastQC counts")

    def upload_picard_metric(self, alignment_id, flowcell_lane_id, aggregation_id, filename, metric_name):

        if not self.picard_metrics:
            self.picard_metrics = self.get_picard_metrics()

        picard_metric  = None
        try:
            picard_metric = open(filename, 'r').read()
        except:
            log.error("Could not read picard metric file %s" % filename)
            return None

        log.debug("Uploading metric contents from: %s" % filename)
        log.debug(picard_metric)

        if not metric_name in self.picard_metrics:
            log.error("Could not find metrics type %s" % metric_name)
            return False

        metric = self.picard_metrics[metric_name]

        if alignment_id:
            object_id = alignment_id
            if not self.alignment_contenttype:
                self.get_alignment_contenttype()
            content_type = self.alignment_contenttype
            lane_info = self.get_flowcell_lane(flowcell_lane_id)

            if not lane_info:
                return False

            label = "FC%s %s %s %s %s" % (lane_info['flowcell_label'],
                lane_info["samplesheet_name"], str(lane_info["lane"]),
                lane_info["barcode_index"], metric_name)
        elif aggregation_id:
            object_id = aggregation_id
            if not self.aggregation_contenttype:
                self.get_aggregation_contenttype()
            content_type = self.aggregation_contenttype
            aggregation_info = self.get_aggregation(aggregation_id)
            log.debug(aggregation_info)
            library_info = self.get_by_full_url(aggregation_info['library'])
            if library_info:
                log.debug(library_info)
            else:
                log.error("Could not fetch %s" % aggregation_info['library'])
                return False
            label = "AGG%d LN%d %s" % (aggregation_id, library_info['number'], metric_name)

        # does this report already exist?
        log.debug("Checking for existing report...")
        existing = self.get_single_result(
            'picard_report/',
            query={
                "object_id": object_id,
                "content_type": content_type['id'],
                "metric": metric['id'],
                "label": label,
            })

        if existing and 'raw_data' in existing and existing['raw_data'] == picard_metric:
            log.info("Picard report is the same, not uploading")
            return

        upload = dict()

        upload['metrics'] = [metric['url']]
        upload['raw_data'] = picard_metric
        upload['content_type'] = content_type["url"]
        upload['object_id'] = object_id
        upload['label'] = label

        if existing is not None:
            result = self.patch(url=existing['url'], json=upload)
        else:
            log.info("Uploading new picard report %s" % upload['label'])
            result = self.post("picard_report/", json=upload)

        if not result:
            log.error("Could not upload new Picard report %s" % filename)
        else:
            log.debug(result)

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

    for fastqc_file in poptions.fastqc_files:
        uploader.upload_fastqc(poptions.flowcell_lane_id, fastqc_file)

    if poptions.fastqc_files and poptions.fastqc_counts:
        uploader.upload_fastqc_counts(poptions.alignment_id)

    if poptions.inserts_file:
        uploader.upload_picard_metric(poptions.alignment_id, poptions.flowcell_lane_id, poptions.aggregation_id, poptions.inserts_file, "CollectInsertSizeMetrics")

    if poptions.dups_file:
        uploader.upload_picard_metric(poptions.alignment_id, poptions.flowcell_lane_id, poptions.aggregation_id, poptions.dups_file, "MarkDuplicates")

    if poptions.spot_file or poptions.spot_dup_file:
        uploader.upload_spot(poptions.alignment_id, poptions.spot_file, poptions.spot_dup_file)

    if poptions.counts_file:
        uploader.upload_counts(poptions.alignment_id, poptions.counts_file)

    if poptions.rna_file:
        uploader.upload_rna_metrics(poptions.alignment_id, poptions.rna_file)

    if poptions.barcode_report_file:
        uploader.upload_barcode_report(poptions.barcode_report_file)

    if poptions.alignment_id and poptions.clear_align_stats:
        uploader.clear_alignment_stats(poptions.alignment_id)

    if poptions.alignment_id and (poptions.version_file or poptions.adapter_file or poptions.align_start_time or poptions.align_complete_time):
        uploader.upload_alignment_records(poptions.alignment_id,
            version_file=poptions.version_file, adapter_file=poptions.adapter_file,
            start_time=poptions.align_start_time, complete_time=poptions.align_complete_time)

    if poptions.aggregation_id and poptions.clear_aggregation_stats:
        uploader.clear_aggregation_stats(poptions.aggregation_id)

    if poptions.aggregation_id and poptions.start_aggregation:
        uploader.start_aggregation(poptions.aggregation_id)

    if poptions.aggregation_id and poptions.complete_aggregation:
        uploader.complete_aggregation(poptions.aggregation_id)

    if poptions.flowcell:
        uploader.clear_flowcell_cache(poptions.flowcell)

    if poptions.attach_file:
        uploader.upload_file_attachment(poptions.attach_file, poptions.attach_file_contenttype, poptions.attach_file_objectid,
            file_type=poptions.attach_file_type, file_purpose=poptions.attach_file_purpose, skip_md5_check=poptions.skip_md5_check)

    if poptions.attach_directory:
        uploader.upload_directory_attachment(poptions.attach_directory, poptions.attach_file_contenttype, poptions.attach_file_objectid,
            file_purpose=poptions.attach_file_purpose)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
