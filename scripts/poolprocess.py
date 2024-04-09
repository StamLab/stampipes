#import csv
import argparse
import functools
import json
import logging
import os
import re
import sys
import textwrap
from collections import OrderedDict, defaultdict

import requests

try:
    from concurrent.futures import ThreadPoolExecutor
except ImportError:
    from futures import ThreadPoolExecutor

# Globals for storing our mapping (saves LIMS hits)
POOL_KEY_TO_LIB_IDS = defaultdict(list)  # {(pool_id, lane_number): [lib_id]}
LIB_ID_TO_LANE_IDS = defaultdict(list)   # {lib_id:                 [lane_ids]}
LANE_ID_TO_ALN_IDS = defaultdict(list)   # {lane_id:                [aln_ids]}
LANES_WITH_DIRECT_POOL = {}


LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Keys to copy directly from the api/locus/ endpoint
LOCUS_KEYS = ["genome_label", "chromosome_name", "genes", "name",
              "genomic_feature", "genomic_coordinate_genome_label",
              "genomic_coordinate_chromosome_name", "genomic_coordinate_start",
              "genomic_coordinate_end"]

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')

SCRIPT_OPTIONS = {
    "quiet": False,
    "debug": False,
    "base_api_url": None,
    "token": None,
    "flowcell": None,
    "tag": None,
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_queue": "hpcz-2",
    "qsub_prefix": ".proc",
    "dry_run": False,
    "no_mask": False,
    "redo_completed": False,
    "qsub_priority": -50,
    "auto_aggregate": False,
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

    #parser.add_argument("--alignment", dest="align_ids", type=int, action="append",
    #    help="Run for this particular alignment.")
    parser.add_argument("--flowcell", dest="flowcell_label",
                        help="Run for this particular flowcell label.")
    #parser.add_argument("--pool", dest="pool",
    #                    help="Run for this particular pool.")
    #parser.add_argument("--tag", dest="tag",
    #    help="Run for alignments tagged here.")
    #parser.add_argument("--project", dest="project",
    #    help="Run for alignments in this project.")

    parser.add_argument("--script_template", dest="script_template",
                        help="The script template to use.")
    parser.add_argument("--qsub_priority", dest="qsub_priority", type=int,
                        help="The priority to give scripts we are submitting.")

    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="Append commands to run this alignment to this file.")
    parser.add_argument("-b", "--sample-script-basename", dest="sample_script_basename",
                        help="Name of the script that goes after the sample name.")
    parser.add_argument("--qsub-prefix", dest="qsub_prefix",
                        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.")
    parser.add_argument("--qsub-queue", dest="qsub_queue",
                        help="Name of the SLURM partition")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
                        help="Take no action, only print messages.")
    parser.add_argument("--no-mask", dest="no_mask", action="store_true",
                        help="Don't use any barcode mask.")
    parser.add_argument("--redo_completed", dest="redo_completed", action="store_true",
                        help="Redo alignments marked as completed.")
    #parser.add_argument("--auto_aggregate", dest="auto_aggregate", help="Script created will also run auto-aggregations after alignments finished.",
        #action="store_true")
    parser.add_argument("--align_base_dir", dest="align_base_dir",
                        help="Create the alignment directory in this directory")

    parser.add_argument("--listout", dest="simple_output", action="store_true",
                        help="Write only a list of alignments to run, rather than a script to submit them")

    parser.set_defaults(**SCRIPT_OPTIONS)
    parser.set_defaults(quiet=False, debug=False)

    return parser


class ProcessSetUp(object):

    def __init__(self, args, api_url, token):

        self.token = token
        self.api_url = api_url
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.dry_run = args.dry_run
        self.no_mask = args.no_mask
        self.redo_completed = args.redo_completed
        self.script_template = args.script_template
        self.qsub_priority = args.qsub_priority
        self.qsub_queue = args.qsub_queue
        #self.auto_aggregate = args.auto_aggregate
        self.align_base_dir = args.align_base_dir

        self.simple_output = args.simple_output

        self.session = requests.Session()
        self.session.headers.update({'Authorization': "Token %s" % self.token})

        self.pool = ThreadPoolExecutor(max_workers=10)

    @functools.lru_cache(maxsize=None)
    def api_single_result(self, url_addition=None, url=None):

        if url_addition:
           url = "%s/%s" % (self.api_url, url_addition)

        request = self.session.get(url)

        if request.ok:
            logging.debug(request.json())
            return request.json()
        else:
            logging.error("Could not get data from %s" % url)
            logging.error(request)
            return None

    @functools.lru_cache(maxsize=None)
    def api_list_result(self, url_addition=None, url=None):

        more = True
        results = []

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        while more:

            logging.debug("Fetching more results for query %s" % url)

            request = self.session.get(url)

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

    def get_align_process_info(self, alignment_id):

        process_info = self.api_single_result("flowcell_lane_alignment/%d/processing_information/" % alignment_id)

        if not process_info:
            logging.critical("Could not find processing info for alignment %d\n" % alignment_id)
            logging.critical(process_info)
            sys.exit(1)

        return process_info

    def get_process_template(self, align_id, process_template_id):

        if not process_template_id:
            logging.critical("No process template for alignment %d\n" % align_id)
            return None

        info = self.api_single_result("process_template/%d/" % (process_template_id))

        if not info:
            logging.critical("Could not find processing template for ID %d\n" % process_template_id)
            sys.exit(1)

        return info

    # Run alignment setup in parallel
    def setup_alignments(self, align_ids, parallel=True):
        all_okay = True
        logging.info("Setting up alignments: %s", align_ids)
        if parallel:
            for id, error in self.pool.map(self.setup_alignment, align_ids):
                if error:
                    logging.error("ALN%d result received, error: %s" % (id, error))
                    all_okay = False
                else:
                    logging.debug("ALN%d result received, OK" % id)
            if not all_okay:
                #logging.critical("Errors during setup, exiting")
                logging.error("Errors during setup, but continuing with other alignments")
        # Sequential version, helpful for debugging
        else:
            for aln_id in align_ids:
                self.setup_alignment(aln_id)

    def setup_alignment(self, align_id):

        try:
            processing_info = self.get_align_process_info(align_id)
            alignment = self.api_single_result("flowcell_lane_alignment/%d/" % (align_id))

            if self.redo_completed or not alignment['complete_time']:
                self.create_script(processing_info, alignment["id"])
                return (align_id, None)
            else:
                logging.info("Skipping completed alignment %d" % align_id)
                return (align_id, None)
        except Exception as e:
            logging.exception("Could not set up alignment %s}: (%s)" % (align_id, e))
            return (align_id, e)

    def get_lane_file(self, lane_id, purpose):
        candidates = self.api_list_result("file/?content_type=40&purpose__slug=%s&object_id=%d" % (purpose, lane_id))

        if not candidates:
            return None
        if len(candidates) > 1:
            return None

        return candidates[0]

    def setup_tag(self, tag_slug):

        align_tags = self.api_list_result("tagged_object/?content_type=47&tag__slug=%s" % tag_slug)

        self.setup_alignments([align_tag["object_id"] for align_tag in align_tags])

    def setup_project(self, project_id):
        logging.info("Setting up project #%s" % project_id)
        alignments = self.api_list_result("flowcell_lane_alignment/?lane__sample__project=%s" % project_id)
        self.setup_alignments([alignment["id"] for alignment in alignments])

    def setup_flowcell(self, flowcell_label):
        logging.info("Setting up flowcell for %s" % flowcell_label)
        align_ids = self.get_alignment_ids(flowcell_label)

        logging.debug("align ids: %s", align_ids)
        #alignments = self.api_list_result("flowcell_lane_alignment/?lane__flowcell__label=%s&page_size=1000" % flowcell_label)
        # Disable parallelism so that caching works
        self.setup_alignments(align_ids, parallel=False)
        self.add_stats_upload(flowcell_label)

    def get_alignment_ids(self, flowcell_label: str) -> [int]:
        """
        For each librarypool/lane combination on the flowcell:
        Pick one representative alignment (the one with the lowest ID)
        return the IDs of those flowcells
        """

        def extract_id_from_url(url):
            if url is None:
                return None
            return int(re.findall(r'\d+', url)[-1])

        # Storage for the 3 layers of mapping between alignments and pools
        global POOL_KEY_TO_LIB_IDS
        global LIB_ID_TO_LANE_IDS
        global LANE_ID_TO_ALN_IDS
        global LANES_WITH_DIRECT_POOL

        POOL_KEY_TO_LIB_IDS = defaultdict(list)  # {(pool_id, lane_number): [lib_id]}
        LIB_ID_TO_LANE_IDS = defaultdict(list)   # {lib_id:                 [lane_ids]}
        LANE_ID_TO_ALN_IDS = defaultdict(list)   # {lane_id:                [aln_ids]}

        library_info = set()
        for lane in self.api_list_result("flowcell_lane/?flowcell__label=%s&page_size=1000" % flowcell_label):
            lib_url = lane['library']
            lane_lane = lane['lane']
            if lib_url is not None:
                library_info.add((lib_url, lane_lane))
                lib_id = extract_id_from_url(lib_url)
                LIB_ID_TO_LANE_IDS[lib_id].append(lane['id'])
            else:
                # HACKS BELOW
                # Get pool manually
                pool_url = lane['library_pool']
                pool_id = extract_id_from_url(pool_url)
                LANES_WITH_DIRECT_POOL[lane['id']] = pool_id
                pool_key = (pool_id, lane_lane)
                pool_number = int(lane['library_pool__number'])

                # Get Library info
                lp_info = self.api_single_result(url=pool_url)
                sl_info = self.api_single_result(url=lp_info['sublibrary'])
                cl_info = self.api_single_result(url=sl_info['cell_library'])
                lib_ids = [extract_id_from_url(lib_url) for lib_url in cl_info["libraries"]]

                for lib_url in cl_info["libraries"]:
                    library_info.add((lib_url, lane_lane))

                for lib_id in lib_ids:
                    POOL_KEY_TO_LIB_IDS[pool_key].append(lib_id)
                    LIB_ID_TO_LANE_IDS[lib_id].append(lane['id'])


        # Set of poolnums + lane
        pool_info = set()
        for info in library_info:
            lib_info = self.api_single_result(url=info[0])
            for pool in lib_info['librarypools']:
                pool_info.add((pool["number"], info[1]))
                key = (pool["id"], info[1])
                POOL_KEY_TO_LIB_IDS[key].append(lib_info['id'])

        all_alignments = self.api_list_result("flowcell_lane_alignment/?lane__flowcell__label=%s&page_size=1000" % flowcell_label)

        direct_alns = set()
        for aln in all_alignments:
            lane_id = extract_id_from_url(aln['lane'])
            LANE_ID_TO_ALN_IDS[lane_id].append(aln['id'])
            if lane_id in LANES_WITH_DIRECT_POOL:
                direct_alns.add(aln['id'])


        # Find the minimum alignment ID for each pool/lane combination
        lowest_aln_for_pool = {pool_key: None for pool_key in POOL_KEY_TO_LIB_IDS.keys()}
        for (pool_key, lib_ids) in POOL_KEY_TO_LIB_IDS.items():
            for lib_id in lib_ids:
                for lane_id in LIB_ID_TO_LANE_IDS[lib_id]:
                    for aln_id in LANE_ID_TO_ALN_IDS[lane_id]:
                        cur_aln = lowest_aln_for_pool[pool_key]
                        logging.debug("%s, %d, %d, %d < %s?",
                                      pool_key, lib_id, lane_id, aln_id, cur_aln)
                        if cur_aln is None or cur_aln > aln_id:
                            lowest_aln_for_pool[pool_key] = aln_id

        align_ids = set(lowest_aln_for_pool.values()).union( direct_alns )
        logging.debug("POOL_KEY_TO_LIB_IDS %s", POOL_KEY_TO_LIB_IDS)
        logging.debug("LIB_ID_TO_LANE_IDS %s", LIB_ID_TO_LANE_IDS)
        logging.debug("LANE_ID_TO_ALN_IDS %s", LANE_ID_TO_ALN_IDS)
        logging.debug("ALN IDS %s", align_ids)
        return list(align_ids)




    #def auto_aggregation_script(self,flowcell_label,alignments):
    #    aaname_sentinel = "auto_agg_sentinel.%s" % (flowcell_label)

    #    if not self.outfile:
    #        logging.debug("Writing script to stdout")
    #        outfile = sys.stdout
    #    else:
    #        logging.debug("Logging script to %s" % self.outfile)
    #        outfile = open(self.outfile, 'a')

    #    contents = textwrap.dedent("""\
    #               cd "$FLOWCELLS"/FC{label}_*
    #               sentinel_dependencies=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterany/g')
    #               sbatch --export=ALL -J {job_name} -o {job_name}.o%A -e {job_name}.e%A --partition={queue} --cpus-per-task=1 --ntasks=1 $sentinel_dependencies --mem-per-cpu=1000 --parsable --oversubscribe <<__AUTOAGG1__
    #               #!/bin/bash
    #               rm -f run_aggregations.bash
    #               python $STAMPIPES/scripts/aggregateprocess.py --flowcell {label} --outfile run_aggregations.bash --qsub-queue {qqueue}
    #               bash run_aggregations.bash
    #               __AUTOAGG1__
    #               """.format(label=flowcell_label,
    #                        job_name=aaname_sentinel,
    #                        queue=self.qsub_queue,
    #                        qqueue=self.qsub_queue))

    #    outfile.write(contents)
    #    outfile.close()

    def add_script(self, align_id, processing_info, script_file, sample_name):

        ram_megabytes = 4000

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, 'a')

        if self.simple_output:
            outfile.write(script_file + "\n")
        else:
            outfile.write("cd %s && " % os.path.dirname(script_file))
            fullname = "%s%s-%s-ALIGN#%d" % (self.qsub_prefix,sample_name,processing_info['flowcell']['label'],align_id)
            outfile.write("jobid=$(sbatch --export=ALL -J %s -o %s.o%%A -e %s.e%%A --partition=%s --cpus-per-task=1 --ntasks=1 --mem-per-cpu=%d --parsable --oversubscribe <<__ALIGNPROC__\n#!/bin/bash\nbash %s\n__ALIGNPROC__\n)\nPROCESSING=\"$PROCESSING,$jobid\"\n\n" % (fullname, fullname, fullname, self.qsub_queue, ram_megabytes, script_file))
        outfile.close()


    def add_stats_upload(self, flowcell_label):
        job_name = ".upload-altcode-%s" % flowcell_label
        template = textwrap.dedent(
            """\
            cd "$FLOWCELLS"/FC{label}_*
            sentinel_dependencies=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterany/g')
            sbatch --export=ALL -J {job_name} -o {job_name}.o%A -e {job_name}.e%A --partition={queue} --cpus-per-task=1 --ntasks=1 $sentinel_dependencies --mem-per-cpu=1000 --parsable --oversubscribe <<__UPLOAD_POOL_DATA__
            #!/bin/bash
            python $STAMPIPES/scripts/altcode/upload_stats.py "$PWD"
            __UPLOAD_POOL_DATA__""")
        content = template.format(
            label=flowcell_label,
            job_name=job_name,
            queue=self.qsub_queue,
        )

        with open(self.outfile, 'a') as outfile:
               outfile.write(content)

    def get_script_template(self, process_template):

        if self.script_template:
            script_path = self.script_template
        else:
            script_path = os.path.expandvars(process_template["process_version"]["script_location"])
        return open(script_path, 'r').read()

    def create_script(self, processing_info, align_id):
        logging.debug("Creating script for ALN%d", align_id)
        assert len(processing_info["libraries"]) == 1
        lane = processing_info["libraries"][0]
        alignment = [a for a in lane["alignments"] if a["id"] == align_id][0]

        if not "process_template" in alignment:
            logging.error("Alignment %d has no process template" % align_id)
            return False

        process_template = self.get_process_template(align_id, alignment["process_template"])

        if not process_template:
            return False

        flowcell_directory = processing_info['flowcell']['directory']

        share_dir = lane.get("project_share_directory")
        if share_dir:
            flowcell_directory = os.path.join(share_dir, "alignments")
        if not flowcell_directory:
            logging.error("Alignment %d has no flowcell directory for flowcell %s" % (align_id, processing_info['flowcell']['label']))
            return False

        if lane.get('library'):
            lib_info_response = self.api_single_result("library/?number=%d" % int(lane["library"]))["results"]
            assert len(lib_info_response) == 1
            lib_info = lib_info_response[0]
            logging.debug("lib info is %s", lib_info)
            pool_name = lib_info["librarypools"][0]["object_name"]
            logging.debug("pool is %s", pool_name)
        else:
            pool_name = lane['samplesheet_name']

        fastq_directory = os.path.join(flowcell_directory, "Project_%s" % lane['project'], "LibraryPool_%s" % pool_name)

        # Reset the alignment's sample name if we decied not to use the barcode index mask
        if self.no_mask:
            alignment['sample_name'] = "%s_%s_L00%d" % (lane['samplesheet_name'], lane['barcode_index'], lane['lane'])

        align_dir = "align_%d_%s_%s" % (alignment['id'], alignment['genome_index'], alignment['aligner'])
        if alignment['aligner_version']:
            align_dir = "%s-%s" % (align_dir, alignment['aligner_version'])

        script_directory = os.path.join(fastq_directory, align_dir)
        if self.align_base_dir:
            script_directory = os.path.join(self.align_base_dir, align_dir)

        r1_fastq = self.get_lane_file(lane["id"], "r1-fastq")

        if not r1_fastq:
            logging.error("Missing r1-fastq for lane %d (alignment %d) - check dir %s" % (lane["id"], alignment["id"], fastq_directory))
            return False

        if processing_info['flowcell']['paired_end']:
            r2_fastq = self.get_lane_file(lane["id"], "r2-fastq")
            if not r2_fastq:
                logging.error("Missing r2-fastq for lane %d (alignment %d)" % (lane["id"], alignment["id"]))
                return False

        script_file = os.path.join( script_directory, "%s-%s" % (alignment['sample_name'], self.qsub_scriptname) )
        logging.info("Will write to %s" % script_file)


        # Set up & add environment variables
        env_vars = OrderedDict()

        env_vars["SAMPLE_NAME"] = alignment['sample_name']
        env_vars["BWAINDEX"]    = alignment['genome_index_location']
        env_vars["GENOME"]      = alignment['genome_index']
        env_vars["ASSAY"]       = lane['assay']
        env_vars["READLENGTH"]  = processing_info['flowcell']['read_length']
        try:
            env_vars["LIBRARY_KIT"] = '"' + processing_info['libraries'][0]['library_kit_method'] + '"'
        except:
            env_vars["LIBRARY_KIT"] = None

        if processing_info['flowcell']['paired_end']:
            env_vars["PAIRED"] = "True"
        else:
            env_vars["PAIRED"] = None

        env_vars["FLOWCELL_LANE_ID"] = lane['id']
        env_vars["ALIGNMENT_ID"]     = alignment['id']
        env_vars["ALIGN_DIR"]        = os.path.join(fastq_directory, align_dir)
        env_vars["R1_FASTQ"]         = r1_fastq["path"]

        if processing_info['flowcell']['paired_end']:
            env_vars["R2_FASTQ"] = r2_fastq["path"]

        env_vars["FASTQ_DIR"] = fastq_directory
        env_vars["FLOWCELL"]  = processing_info['flowcell']['label']

        if "barcode1" in lane and lane["barcode1"]:
            p7_adapter = lane['barcode1']['adapter7']
            p5_adapter = lane['barcode1']['adapter5']
            if "barcode2" in lane and lane['barcode2']:
                # Override the "default" end adapter from barcode1
                p5_adapter = lane['barcode2']['adapter5_reverse_complement']

            if not p7_adapter or not p5_adapter:
                logging.warn("Alignment %d missing adapters, some processes might not work" % alignment['id'])

            env_vars["ADAPTER_P7"] = p7_adapter
            env_vars["ADAPTER_P5"] = p5_adapter

            # Process with UMI if the barcode has one and this is a dual index
            # flowcell
            if lane['barcode1']['umi'] and processing_info['flowcell']['dual_index']:
                env_vars["UMI"] = "True"
            else:
                env_vars["UMI"] = None
            env_vars["UMI_METHOD"] = lane['barcode1']['umi_method']

        # Set process template env var overrides
        if 'process_variables' in process_template and process_template['process_variables']:
            try:
                process_template_variables = json.loads(process_template['process_variables'],
                                                        object_pairs_hook=OrderedDict)
                for var, value in process_template_variables.items():
                    env_vars[var] = value
            except ValueError as e:
                logging.error("Could not parse process variables for align %d (template %d): '%s'" %
                              (
                                  alignment['id'],
                                  process_template['id'],
                                  process_template['process_variables']
                              ))
                return False

        if self.dry_run:
            logging.info("Dry run, would have created: %s" % script_file)
            logging.debug(env_vars)
            self.create_sample_config(processing_info, alignment, script_directory, pool_name)
            return True

        if not os.path.exists(script_directory):
            logging.info("Creating directory %s" % script_directory)
            os.makedirs(script_directory)

        # Append to master script
        self.add_script(align_id, processing_info, script_file, alignment['sample_name'])

        # Write file
        outfile = open(script_file, 'w')
        outfile.write("set -e -o pipefail\n")

        # Set env vars
        for var, value in env_vars.items():
            if value is not None:
                outfile.write("export %s=%s\n" % (var, value))
            else:
                outfile.write("unset %s\n" % var)

        outfile.write("\n")
        outfile.write("export QUEUE=%s\n" % (self.qsub_queue))
        outfile.write("\n")
        outfile.write(self.get_script_template(process_template))
        outfile.close()

        # Create the config file as well
        self.create_sample_config(processing_info, alignment, script_directory, pool_name)

    def create_sample_config(self, processing_info, alignment, script_directory, pool_name):
        alignment_id = int(alignment["id"])
        logging.debug("Creating sample config for ALN%d", alignment_id)

        def get_libraries_in_pool(alignment_id):

            # Get all lane ids
            # Go up to the pool then down to the lanes
            # Note: This is inefficient but probably doesnt matter in practice
            lanes = []
            lanes_with_align = set()
            for (lane_id, aln_ids) in LANE_ID_TO_ALN_IDS.items():
                if alignment_id in aln_ids:
                    lanes_with_align.add(lane_id)
            assert len(lanes_with_align) == 1, "Alignment must have exactly 1 lane"
            align_lane_id = lanes_with_align.pop()

            libs_with_align = set()
            for (lib_id, lane_ids) in LIB_ID_TO_LANE_IDS.items():
                if align_lane_id in lane_ids:
                    libs_with_align.add(lib_id)
            #assert len(libs_with_align) == 1, "Lane must have exactly 1 library"
            align_lib_id = libs_with_align.pop()

            pools_with_align = set()
            for (pool_key, lib_ids) in POOL_KEY_TO_LIB_IDS.items():
                if align_lib_id in lib_ids:
                    pools_with_align.add(pool_key)
            # TODO: This is broken because the pool can be in more than one lane!!!
            #assert len(pools_with_align) == 1, "Lib must have exactly one pool"
            align_poolkey = pools_with_align.pop()
            logging.debug("Alignment ALN%d - poolkey %s", alignment_id, align_poolkey)

            library_ids = set(POOL_KEY_TO_LIB_IDS[align_poolkey])
            logging.debug("Lib IDs in poolkey %s: %s", align_poolkey, library_ids)
            return library_ids

        lib_ids = get_libraries_in_pool(alignment_id)

        def build_library_info(lib_id, flowcell_label):
            errors = []
            def add_error(fmt, *args):
                err_msg = fmt % args
                errors.append(err_msg)
                logging.error(fmt, *args)

            lib_info = self.api_single_result("library/%d/" % lib_id)
            barcode = ""
            bc1 = lib_info["barcode1__sequence"]
            bc2 = lib_info["barcode2__sequence"]
            if bc1 is not None:
                barcode += bc1
            if bc2 is not None:
                barcode += bc2

            sample_info = self.api_single_result(url=lib_info["sample"])
            project_info = self.api_single_result(url=sample_info["project"]) if sample_info.get("project") else {"name": None}

            taggedobject_infos = self.api_list_result("tagged_object/?object_id=%d&content_type=%d"
                                             % (lib_info["id"], lib_info["object_content_type"]))
            cycle = None
            for taggedobject_info in taggedobject_infos:
                # TODO: It may be better to check membership in the Insights tag
                if taggedobject_info["tag_slug"].startswith("megamap-run-mmap") or taggedobject_info["tag_slug"].startswith("epicapdev-run-ecd"):
                    if cycle is None:
                        tag_slug = str(taggedobject_info["tag_slug"])
                        match = re.search(r"\d+$", tag_slug)
                        if match:
                            cycle = int(match.group())
                        else:
                            add_error("problem tag slug is '%s'", tag_slug)
                    else:
                        add_error("Multiple tags for LN%d", lib_info["number"])

            def build_effector_info(effectortopool):
                eff = effectortopool["assemble_effector"]
                talen_number = eff["talen_sheet_number"]
                talen = "TL%s" % talen_number if talen_number else None

                return {
                    "chromosome": eff["chromosome"],
                    "start": eff["start"],
                    "end": eff["end"],
                    "strand": eff["strand"],
                    "working_name": eff["working_name"],
                    "talen": talen,

                    "n_terminus": {
                        "name": eff["effector__n_terminus__name"],
                        "nucleotide": eff["effector__n_terminus__nucleotide"],
                        "functional_domain": eff["effector__n_functional_domain__name"],
                    },
                    "c_terminus": {
                        "name": eff["effector__c_terminus__name"],
                        "nucleotide": eff["effector__c_terminus__nucleotide"],
                        "functional_domain": eff["effector__c_functional_domain__name"],
                    },
                    "concentration": eff["concentration"],
                    "ratio_260to280": eff["ratio_260to280"],
                    "ivt_mrna_concentration": eff["ivt_mrna_concentration"],
                    "repeat_array__target_recognition_sequence": "CTCTTTCACAGCTCGCG",
                    "target_recognition_sequence": "TCTCTTTCACAGCTCGCGT",
                    "wells": [
                        {
                            "plate_name": well["plate__name"],
                            "plate_id": well["plate_id"],
                            "well": well["label"],
                        }
                        for well in eff["plate_wells"]
                    ]
                }

            pool_info = []
            # Dummy info in case we can't get it from LIMS
            tc_info = {"notes": "", "number": 0, "sample_taxonomy__name": ""}
            try:
                tc_info = self.api_single_result(url=sample_info["tissue_culture"])
                for effector_pool in tc_info["effector_pools"]:
                    effector_pool_info = self.api_single_result(url=effector_pool["url"])
                    loci_info = []
                    if effector_pool_info.get("loci", False):
                        for locus_url in effector_pool_info["loci"]:
                            locus_info = self.api_single_result(url=locus_url)
                            locus_dict = {
                                "label": locus_info.get("object_label"),
                            }
                            for key in LOCUS_KEYS:
                                locus_dict[key] = locus_info.get(key, None)
                            loci_info.append(locus_dict)

                    pool_info.append({
                        "effector_pool": effector_pool_info["object_name"],
                        "name": effector_pool_info["name"],
                        "purpose": effector_pool_info["purpose__name"],
                        "loci": loci_info,
                        "effectors": [
                            build_effector_info(efftopool)
                            for efftopool in effector_pool_info["effectortopool_set"]
                        ],
                    })
            except:
                add_error("Could not get effector information for sample DS%s", sample_info['number'])

            def extract_lenti_from_tc_notes(notes):
                def match_notes(regex):
                    match = re.search(regex, notes, re.MULTILINE | re.IGNORECASE)
                    if match is None:
                        return None
                    return match.group(1)
                # Example notes field below:
                # Talen Number: TL120935
                # Original TALE name: IL2RA-TL52068-Z
                # LentiTALE: Lenti-KRAB
                # MOI Estimate: 1.4
                # Virus volume: 100
                # Lenti-X Content: 15%
                talen_number = match_notes(r"Talen Number: (TL\d+)\s*")
                talen_name = match_notes(r"Original TALE name: (.+?)\s*$")
                lentitale = match_notes(r"LentiTALE: (.+?)\s*$")
                moi_estimate = match_notes(r"MOI Estimate: (.+?)\s*$")
                virus_volume = match_notes(r"Virus volume: (\d+)\s*$")
                lenti_x_content = match_notes(r"Lenti-X Content: (.+?)\s*$")
                effector_assembly_qc = match_notes(r"Effector Assembly QC: (.+?)\s*$")

                return {
                    "talen_number": talen_number,
                    "talen_name": talen_name,
                    "lentitale": lentitale,
                    "moi_estimate": moi_estimate,
                    "virus_volume": virus_volume,
                    "lenti_x_content": lenti_x_content,
                    "effector_assembly_qc": effector_assembly_qc,
                }

            def sample_plate_wells(sample_info) -> "[dict]":
                def info_to_data(well_info):
                    match = re.match(r"(.*) ([A-Z0-9]{2})", well_info["object_name"])
                    well_data = {
                        "plate_name": match.group(1),
                        "well_label": well_info["label"],
                        "plate_id": well_info["plate_id"],
                        "object_label": well_info["content_object_label"],
                    }
                    if well_info["parent"]:
                        parent_info = self.api_single_result(url=well_info["parent"])
                        well_data["well_parent"] = info_to_data(parent_info)
                    return well_data
                wells = []
                for well in sample_info["plate_wells"]:
                    well_info = self.api_single_result("plate_well/%d/" % well["id"])
                    well_data = info_to_data(well_info)
                    wells.append(well_data)
                return wells

            def reverse_complement(bc: "Optional[str]") -> "Optional[str]":
                if bc is None:
                    return None
                lookup = {"A":"T", "T":"A", "C":"G", "G":"C"}
                return "".join(reversed([lookup[c] for c in bc]))

            lenti_from_tc = extract_lenti_from_tc_notes(tc_info["notes"])
            lib_plate_wells = sample_plate_wells(lib_info)
            deep_info = {
                "barcode": barcode,
                "barcode1": bc1,
                "barcode2": bc2,
                "library": "LN%d" % lib_info["number"],
                "sublibrary": lib_info["sub_library"],
                "sample": "DS%d" % lib_info["sample_number"],
                "tc": "TC%d" % tc_info["number"],
                "tc_notes": tc_info["notes"],
                "lentitale_from_tc_notes": lenti_from_tc,
                "cell_type": tc_info["sample_taxonomy__name"],
                "library_plate_wells": lib_plate_wells,
                "sample_plate_wells": sample_plate_wells(sample_info),
                "project": project_info["name"],
                "flowcell": flowcell_label,
                "cycle": cycle,
                "effector_pools": pool_info,
                "pool": pool_name,
            }

            #####################################
            # Refine our data to match the spec #
            #####################################
            seq_well_label = None
            seq_well_plate = None
            sample_well_label = None
            sample_well_plate = None
            tc_well_label = None
            tc_well_plate = None
            try:
                seq_well_label    = lib_plate_wells[0]["well_label"]
                seq_well_plate    = "PL%d" % lib_plate_wells[0]["plate_id"]
                sample_well_label = lib_plate_wells[0]["well_parent"]["well_label"]
                sample_well_plate = "PL%d" % lib_plate_wells[0]["well_parent"]["plate_id"]
                tc_well_label     = lib_plate_wells[0]["well_parent"]["well_parent"]["well_label"]
                tc_well_plate     = "PL%d" % lib_plate_wells[0]["well_parent"]["well_parent"]["plate_id"]
            except Exception as e:
                add_error("Could not find well info in %s", lib_plate_wells)

            def sort_talens(tls):
                """ Sort talens by number """
                def get_num(tl):
                    match = re.search(r"TL(\d+)", tl)
                    if match:
                        return int(match.group(1))
                    else:
                        logging.warning("Weird talen: '%s'" % tl)
                        return 0
                return sorted(tls, key=get_num)

            if pool_info:
                #talen_name = None #TODO
                talen_names = []
                for pool in deep_info["effector_pools"]:
                    for effector in pool.get("effectors", []):
                        if effector["talen"]:
                            talen_names.append(effector["talen"])

                talen_name = ",".join(sort_talens(talen_names))
                orig_talen_name = talen_name

            else:
                talens_str = lenti_from_tc["talen_name"]
                orig_talen_name = lenti_from_tc["talen_name"]
                # Sort to normalize for string comparison
                if talens_str is not None:
                    talens = sort_talens([t.strip() for t in talens_str.split(",")])
                    talen_name = ",".join(talens)
                else:
                    talen_name = None

            lenti_qc_passed = lenti_from_tc["effector_assembly_qc"] is None

            if sample_info["time_point_unit"] == 5:
                # harvest timepoint is in days
                harvest_timepoint = float(sample_info["time_point"])
            else:
                add_error("Sample timepoint unit unknown")
                harvest_timepoint = None

            # def parse_talen_names_from_tc_notes(notes):
            #     def match_notes(regex):
            #         match = re.search(regex, notes, re.MULTILINE | re.IGNORECASE)
            #         if match is None:
            #             return None
            #         return match.group(1)
            #     talen_new = match_notes(r"Talen Number:\s*(.+?)\s*$")
            #     talen_orig = match_notes(r"Original TALE name:\s*(.+?)\s*$")
            #     if talen_new is not None and talen_orig is not None:
            #         # Sort to make matching easier
            #         split_new = [s.strip() for s in talen_new.split(",")]
            #         split_orig = [s.strip() for s in talen_new.split(",")]
            #         if len(split_new) != len(split_orig):
            #             log.warning("Length of new and old talens differ")
            #             break
            #         sorted_pairs = sorted(zip(split_new, split_orig))
            #         talen_new = [p[0] for p in sorted_pairs]
            #         talen_orig = [p[1] for p in sorted_pairs]
            #     return (talen_orig, talen_new)

            # (talen_orig, talen_new) = parse_talen_names_from_tc_notes(tc_info["notes"])

            def get_sbl_and_cl(pool_name):
                (sbl, cl) = (None, None)
                try:
                    m = re.match(r"LP(\d+)", pool_name)
                    if not m:
                        add_error("Pool name '%s' not valid, can't get SBL&CL", pool_name)
                        return (None, None)
                    pool_id = int(m.group(1))
                    pool_info = self.api_list_result("library_pool/?number=%d" % pool_id)[0]
                    if not pool_info.get("sublibrary"):
                        return (None, None)
                    sbl_info = self.api_single_result(url=pool_info.get("sublibrary"))
                    sbl = sbl_info["object_name"]
                    if not sbl_info.get("cell_library"):
                        return (sbl, None)
                    cl_info = self.api_single_result(url=sbl_info.get("cell_library"))
                    cl = cl_info.get("object_name")
                except Exception as e:
                    add_error("Error finding SBL or CL: %s", e)
                return (sbl, cl)
            (sbl_name, cl_name) = get_sbl_and_cl(pool_name)

            info = {
                "sequencing_barcode_well": seq_well_label,
                "sequencing_barcode_plate": seq_well_plate,
                "sample_well": sample_well_label,
                "sample_plate": sample_well_plate,
                "perturbation_well": tc_well_label,
                "perturbation_plate": tc_well_plate,
                "day": harvest_timepoint,
                "talen_name": talen_name,
                "cell_type": tc_info["sample_taxonomy__name"],
                "cycle": cycle,
                "tale_target_name": "TODO",
                "tale_target_master_gene_id": "TODO",
                "effector_purpose": "TODO",
                "cell_library": cl_name,
                "sublibrary": sbl_name,
                "library_pool": pool_name,
                "TC#": "TC%d" % tc_info["number"],
                "DS#": "DS%d" % sample_info["number"],
                "LN#": "LN%d" % lib_info["number"],
                "TL#_original": orig_talen_name,
                "TL#_new": lenti_from_tc["talen_number"],
                "sample_barcode": reverse_complement(bc2),
                "lenti_qc_passed": True,

                "script_errors": errors,
                "additional_information": deep_info,
            }
            return info


        flowcell_label = "FC%s" % processing_info["flowcell"]["label"]
        libraries = []

        for lib_id in lib_ids:
            libraries.append(build_library_info(lib_id, flowcell_label))

        data = {"libraries": libraries}
        if self.dry_run:
            logging.info("dry_run, would have written %s/pool_info.json", script_directory)
            return
        # do stuff
        with open("%s/pool_info.json" % script_directory, "w") as out:
            json.dump(data, out, indent=2, sort_keys=True)
            #writer = csv.DictWriter(out, fieldnames=fieldnames, dialect="excel-tab", restval="")
            #writer.writeheader()
            #for row in rows:
                #writer.writerow(row)


def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=LOG_FORMAT)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
        logging.getLogger("requests").setLevel(logging.WARNING)

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

    process = ProcessSetUp(poptions, api_url, token)

    #process.setup_alignments(poptions.align_ids)

    if poptions.flowcell_label:
        process.setup_flowcell(poptions.flowcell_label)
    else:
        logging.critical("Non-flowcell setup not yet supported")

    #if poptions.tag:
    #    process.setup_tag(poptions.tag)

    #if poptions.project:
    #    process.setup_project(poptions.project)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
