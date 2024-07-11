#!/usr/bin/env python3

import json
import os
import sys
import argparse
import logging
import re
import copy
import requests
import datetime

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.getLogger("requests").setLevel(logging.WARNING)
util_log = logging.getLogger("StamPy.util")


def foldercheck(*args):
    """Checks to see if the folders exist, creates them if they are not."""
    for folder in args:
        if not os.path.isdir(folder):
            try:
                os.mkdir(folder)
                util_log.info("Created folder: %s" % folder)
            except OSError as x:
                util_log.error("ERROR: Could not create directory: %s" % folder)
                util_log.warn(
                    "Please make sure all nonexistant parent directories have been created."
                )
                sys.exit(0)


options = {
    "quiet": False,
    "debug": False,
    "process_config": "processing.json",
    "trackhub_config": None,
    "projectname": None,
    "priority": None,
    "base_api_url": os.getenv("LIMS_API_URL"),
    "token": os.getenv("LIMS_API_TOKEN"),
}


def mysql_clean(input):
    # Mysql names can contain only 0-9, a-z, A-Z, _, or $
    # So we replace all other characters with an underscore,
    # after removing leading/trailing whitespace.
    output = re.sub("[^\w$]", "_", input.strip())
    return output


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
        "-j",
        "--json-config",
        dest="process_config",
        help="The process config to work off of.",
    )
    parser.add_argument(
        "-c",
        "--trackhub-config",
        dest="trackhub_config",
        help="The trackhub config to work off of.",
    )
    parser.add_argument(
        "-p",
        "--priority",
        dest="priority",
        required=True,
        help="The priority of this project",
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
        "-n", "--project-name", dest="projectname", help="Name of project.  Required."
    )
    parser.set_defaults(**options)
    parser.set_defaults(quiet=False, debug=False)
    return parser


############


class MakeBrowserLoad(object):
    genome_organisms = {
        "hg19": "human",
        "hg38": "human",
        "rn5": "rat",
        "mm9": "mouse",
        "mm10": "mouse",
        "TAIR9": "arabidopsis",
        "sacCer2": "sacCer",
        "sacCer3": "sacCer",
        "ce4": "worm",
        "cb3": "worm",
        "K12": "e.coli",
        "NC_000913.2": "e.coli",
        "hera1": "butterfly",
        "hmel1a": "butterfly",
        "panu2a": "baboon",
        "felCat5": "cat",
        "borrBurg": "bacteria",
        "danRer7": "zebrafish",
    }

    rna_strands = ["all", "pos", "neg"]

    def __init__(
        self,
        group_data,
        trackhubconfig,
        priority,
        projectname,
        date,
        base_api_url,
        token,
    ):
        self.win = 75
        self.binI = 20
        self.bigwig = True

        self.data = group_data
        self.priority = priority
        self.projectname = projectname
        self.date = date
        self.load_config(trackhubconfig)
        self.outdir = os.path.join(self.aggregation_link_folder, self.projectname)
        self.base_api_url = base_api_url
        self.token = token

    def load_config(self, trackhubconfig):
        import configparser

        Config = configparser.ConfigParser()
        Config.read(trackhubconfig)
        self.trackhubURL = Config.get("browser", "trackhub_url")
        self.aggregation_link_folder = Config.get("browser", "aggregation_link_folder")

    def load(self):
        # set up folder
        logging.info("Checking for trackhub folder: %s" % self.aggregation_link_folder)
        foldercheck(self.aggregation_link_folder)
        logging.info("Checking for project trackhub folder: %s" % self.outdir)
        foldercheck(self.outdir)

        # prepare tracks for writing
        self.prepare_tracks()
        self.prepare_genomes()

        # write tracks and hub files
        self.create_ras()
        self.create_hubtxt()
        self.create_genometxt()

    def create_hubtxt(self):
        hubfile = os.path.join(self.outdir, "hub.txt")
        logging.info("Creating hub.txt file: %s" % hubfile)
        hub = open(hubfile, "w")
        hub.write("hub %s\n" % self.projectname)
        hub.write("shortLabel %s\n" % self.projectname)
        hub.write("longLabel %s, %s\n" % (self.projectname, self.date))
        hub.write("genomesFile genomes.txt\n")
        hub.write("email anishida@altius.org\n")
        hub.close()

    def create_genometxt(self):
        genomefile = os.path.join(self.outdir, "genomes.txt")
        logging.info("Creating genome.txt file: %s" % genomefile)
        genomes = open(genomefile, "w")
        for key in self.all_tracks:
            genomes.write("\ngenome %s\n" % key)
            genomes.write("trackDb %s/trackDb.%s.txt\n" % (key, self.projectname))
        genomes.close()

    def prepare_tracks(self):
        self.all_tracks = {}

        for agg in self.data:
            # skip aggregations that are not completed
            if agg["needs_reprocessing"] == True:
                continue
            if agg["processing_completed"] == None:
                continue

            tracks = {}
            logging.debug("Preparing tracks for AGG: %s" % agg["id"])
            tracks["agg_id"] = agg["id"]
            tracks["agg_ln"] = agg["library_name"]
            tracks["agg_taxonomy"] = agg["taxonomy_name"]
            tracks["agg_stat"] = ""

            if "stats" in agg and agg["stats"] is not None:
                if "hotspot1-SPOT" in agg["stats"]:
                    tracks["agg_stat"] = agg["stats"]["hotspot1-SPOT"]

            # get genome name (not in json)
            agg_genome_req = requests.get(
                "%s/genome_index/%s" % (self.base_api_url, agg["genome_index_id"]),
                headers={"Authorization": "Token %s" % self.token},
            )
            if agg_genome_req.ok:
                agg_genome_result = agg_genome_req.json()
            # change sequins to normal
            if agg_genome_result["label"] == "GRCh38_no_alts_sequins":
                agg_genome_result["label"] = "GRCh38_no_alts"
            tracks["agg_genome"] = agg_genome_result["label"]

            # output expected is explicit for type of aggregation template used
            # dna
            if (
                agg["aggregation_process_template_id"] == 5
                or agg["aggregation_process_template_id"] == 43
                or agg["aggregation_process_template_id"] == 44
            ):
                if (
                    "normalized-density-bigwig-windowed" in agg["files"]
                    and "density-bigwig-windowed" in agg["files"]
                    and "all-alignments-bam" in agg["files"]
                    and "cutcounts-bw" in agg["files"]
                ):
                    tracks["dnase_normdens"] = agg["files"][
                        "normalized-density-bigwig-windowed"
                    ]
                    tracks["dnase_dens"] = agg["files"]["density-bigwig-windowed"]
                    tracks["dnase_align"] = agg["files"]["all-alignments-bam"]
                    tracks["dnase_cutconts"] = agg["files"]["cutcounts-bw"]
                else:
                    logging.info("Unable to locate AGG files for: %s" % (agg["id"]))
            # rna (processes are seperate for each genome)
            elif (
                agg["aggregation_process_template_id"] == 30
                or agg["aggregation_process_template_id"] == 31
                or agg["aggregation_process_template_id"] == 35
                or agg["aggregation_process_template_id"] == 37
            ):
                if (
                    "all-alignments-bam" in agg["files"]
                    and "neg-coverage-bigwig" in agg["files"]
                    and "pos-coverage-bigwig" in agg["files"]
                ):
                    tracks["rna_align"] = agg["files"]["all-alignments-bam"]
                    tracks["rna_poscov"] = agg["files"]["pos-coverage-bigwig"]
                    tracks["rna_negcov"] = agg["files"]["neg-coverage-bigwig"]
                else:
                    logging.info("Unable to locate AGG files for: %s" % (agg["id"]))
                # coverage across both strands still new, seperate from the rest for now
                if "all-coverage-bigwig" in agg["files"]:
                    tracks["rna_bothcov"] = agg["files"]["all-coverage-bigwig"]
                else:
                    logging.info(
                        "Unable to locate combined stranded AGG files for: %s"
                        % (agg["id"])
                    )
            else:
                logging.info(
                    "Unknown template type, %s, for %s"
                    % (agg["aggregation_process_template_id"], agg["id"])
                )
            if not tracks["agg_genome"] in self.all_tracks:
                self.all_tracks[tracks["agg_genome"]] = []
            self.all_tracks[tracks["agg_genome"]].append(tracks)

    def prepare_genomes(self):
        # change genome names to match UCSC info
        for key in self.all_tracks:
            if key == "mm10-encode3-male":
                self.all_tracks["mm10"] = self.all_tracks.pop("mm10-encode3-male")
            elif key == "GRCh38_no_alts":
                self.all_tracks["hg38"] = self.all_tracks.pop("GRCh38_no_alts")
            elif key == "dm6_minimal":
                self.all_tracks["dm6"] = self.all_tracks.pop("dm6_minimal")

    def create_ras(self):
        for key in self.all_tracks:
            self.create_ra(key)

    def create_ra(self, genome):
        logging.info("Creating RA file for genome, %s" % genome)

        subtracks = self.all_tracks[genome]

        # collect unique path types per genome and samples
        path_names = {}
        for agg in subtracks:
            for info_type in agg:
                if re.match("dnase*", info_type) or re.match("rna*", info_type):
                    path_names[info_type] = 0
        all_samples = set(subtrack["agg_id"] for subtrack in subtracks)

        # write output strings for views and paths
        view_string = ""
        sample_string = ""
        for path in path_names:
            view_string = view_string + " " + path + "=" + path
        for agg in all_samples:
            sample_string = sample_string + " " + str(agg) + "=AG" + str(agg)

        # create genome folder
        foldercheck(os.path.join(self.outdir, genome))

        # create genome RA file
        ra_file = os.path.join(self.outdir, genome, "trackDb.%s.txt" % self.projectname)
        ra = open(ra_file, "w")

        # get sample list up front
        all_samples = set(subtrack["agg_id"] for subtrack in subtracks)

        # write header
        ra.write("track %s\n" % self.projectname)
        ra.write("compositeTrack on\n")
        ra.write("shortLabel %s\n" % self.projectname)
        ra.write("longLabel %s, %s\n" % (self.projectname, self.date))
        ra.write("group %s\n" % self.projectname)
        ra.write("priority %s\n" % self.priority)
        ra.write("subGroup1 view Views%s\n" % view_string)
        ra.write("subGroup2 show _ show=_\n")
        ra.write("dimensions dimensionX=view dimensionY=show\n")
        ra.write("sortOrder view=+\n")
        ra.write("dragAndDrop subTracks\n")
        ra.write("type bed 3 +\n")
        ra.write("noInherit on\n\n")

        # for each path, for each agg
        for path in path_names:
            logging.info("Writing %s to RA for %s", path, genome)

            # hardcoded display settings
            # change alignments to BAM files
            # change auto display settings for normalized densities and RNA pos/neg densities
            file_format = "bigWig"
            visibility = "hide"
            if path == "dnase_align" or path == "rna_align":
                file_format = "bam"
            if (
                path == "dnase_normdens"
                or path == "rna_poscov"
                or path == "rna_negcov"
                or path == "rna_allcov"
            ):
                visibility = "full\n\tviewLimits 0:5\n\tautoScale off\n\tmaxHeightPixels 100:32:16"

            # write path header
            ra.write("\ttrack %s_%s\n" % (self.projectname, path))
            ra.write("\tsubTrack %s\n" % self.projectname)
            ra.write("\tview %s\n" % path)
            ra.write("\tshortLabel %s\n" % path)
            ra.write("\tvisibility %s\n\n" % visibility)

            # write aggs
            for track in subtracks:
                if path in track:
                    friendly_path = re.sub(
                        "/net/seq/data/", self.trackhubURL, track[path]
                    )
                    ra.write(
                        "\t\ttrack %s_%s_%s\n"
                        % (self.projectname, path, track["agg_id"])
                    )
                    ra.write("\t\tbigDataUrl %s\n" % friendly_path)
                    ra.write("\t\tsubTrack %s_%s\n" % (self.projectname, path))
                    ra.write(
                        "\t\tshortLabel AG%s_%s_%s\n"
                        % (track["agg_id"], self.projectname, path)
                    )
                    ra.write("\t\tsubGroups view=%s show=show\n" % path)
                    ra.write(
                        "\t\tlongLabel AG%s, %s, SPOT1 %s, %s, %s\n"
                        % (
                            track["agg_id"],
                            track["agg_ln"],
                            track["agg_stat"],
                            track["agg_taxonomy"],
                            path,
                        )
                    )
                    ra.write("\t\tgroup %s\n" % self.projectname)
                    if file_format == "bam":
                        ra.write("\t\tpairEndsByName .\n")
                    ra.write("\t\ttype %s\n\n" % file_format)

        ra.close()


def main(args=sys.argv):
    parser = parser_setup()
    global poptions
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
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

    data = json.loads(open(poptions.process_config, "r").read())
    dataresults = data["results"]
    dtime = datetime.datetime.now()
    date = "%s-%s-%s" % (dtime.year, dtime.month, dtime.day)

    hubwriter = MakeBrowserLoad(
        dataresults,
        poptions.trackhub_config,
        poptions.priority,
        poptions.projectname,
        date,
        poptions.base_api_url,
        poptions.token,
    )
    hubwriter.load()


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
