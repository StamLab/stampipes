#!/usr/bin/env python3

import argparse
import json
import logging
import os
import re
import sys

import requests

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.getLogger("requests").setLevel(logging.WARNING)

options = {
    "quiet": False,
    "debug": False,
    "process_config": "processing.json",
    "trackhub_config": None,
    "priority": None,
    "api_url": os.getenv("LIMS_API_URL"),
    "api_token": os.getenv("LIMS_API_TOKEN"),
}

util_log = logging.getLogger("StamPy.util")


def foldercheck(*args):
    """Checks to see if the folders exist, creates them if they are not."""
    for folder in args:
        if not os.path.isdir(folder):
            try:
                os.mkdir(folder)
                util_log.info("Created folder: %s" % folder)
            except OSError:
                util_log.error("ERROR: Could not create directory: %s" % folder)
                util_log.warn(
                    "Please make sure all nonexistant parent directories have been created."
                )
                sys.exit(0)


def mysql_clean(input):
    # Mysql names can contain only 0-9, a-z, A-Z, _, or $
    # So we replace all other characters with an underscore,
    # after removing leading/trailing whitespace.
    output = re.sub(r"[^\w$]", "_", input.strip())
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
        help="The priority of this flowcell",
    )
    parser.add_argument(
        "--pre-align-dir",
        dest="pre_align_dir",
        action="store_true",
        help="This flowcell was made before per-alignment directories",
    )
    parser.set_defaults(**options)
    parser.set_defaults(quiet=False, debug=False)
    return parser


class MakeBrowserload(object):
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
        basedir,
        outdir,
        mersize,
        priority,
        paired_end,
        project,
        label,
        date,
    ):
        self.basedir = basedir
        self.flowcell_date = date
        self.outdir = outdir
        self.mersize = mersize
        self.win = 75
        self.binI = 20
        self.priority = priority
        self.paired_end = paired_end
        self.project = project
        self.bigwig = True
        self.projects = project.split(",")
        self.label = label
        self.date = date
        self.project_dirs = {}
        self.data = group_data
        project_dir = ""

        if len(self.projects) == 1 and project_dir:
            self.project_dirs[project] = project_dir
            logging.info("Using project dir: %s" % self.project_dirs[project])
        else:
            for project in self.projects:
                self.project_dirs[project] = os.path.join(
                    self.basedir, "Project_" + project
                )

        self.load_config(trackhubconfig)

    def load_config(self, trackhubconfig):
        import configparser

        config = configparser.ConfigParser()
        config.read(trackhubconfig)
        self.trackhubURL = config.get("browser", "trackhub_url")
        self.flowcell_link_folder = config.get("browser", "flowcell_link_folder")

    def load(self):
        self.basedir_name = os.path.basename(self.basedir)
        foldercheck(self.outdir)

        if False:
            self.main_label = "%s%son%s" % (self.project, self.maintrackname, self.date)
            self.flowcell_name = self.maintrackname
            self.flowcell_date = self.date
        else:
            match = re.search("(FC[A-Z0-9]+)_([0-9]{6})_tag", self.basedir)

            if not match:
                logging.error("Could not figure out a main track name, no flowcell")
                sys.exit(1)

            self.flowcell_name = match.groups()[0]
            if not self.flowcell_date:
                self.flowcell_date = match.groups()[1]

            logging.info("FLOWCELL DATE: %s" % self.flowcell_date)

            self.main_label = "%s%son%s" % (
                self.project,
                self.flowcell_name,
                self.flowcell_date,
            )

        logging.info("Main track name: %s" % self.main_label)

        self.excludes_file = os.path.join(self.outdir, "excludes.%s" % self.main_label)

        if self.flowcell_link_folder:
            logging.debug(
                "link folder: "
                + self.flowcell_link_folder
                + " base folder: "
                + self.basedir_name
            )
            self.link_dir = os.path.join(self.flowcell_link_folder, self.basedir_name)
        else:
            self.link_dir = ""

        self.prepare_tracks()
        logging.info("Main label: %s" % self.main_label)

        # LIMS records early mm10 alignments as 'mm10-encode3-male'
        # change it back to just 'mm10'
        for key in self.subtrack_sets.keys():
            if key == "mm10-encode3-male":
                self.subtrack_sets["mm10"] = self.subtrack_sets.pop("mm10-encode3-male")
            if key == "GRCh38_no_alts":
                self.subtrack_sets["hg38"] = self.subtrack_sets.pop("GRCh38_no_alts")
            if key == "GRCh38_no_alts_sequins":
                self.subtrack_sets["hg38"] = self.subtrack_sets.pop(
                    "GRCh38_no_alts_sequins"
                )

        self.create_ras()
        self.create_hubtxt()
        self.create_genomestxt()
        self.create_htmls()

    # function for creating hub.txt
    def create_hubtxt(self):
        hubfile = os.path.join(self.outdir, "hub.txt")
        logging.info("Creating hub.txt file: %s" % hubfile)
        hub = open(hubfile, "w")
        hub.write("hub %s\n" % self.flowcell_name)
        hub.write("shortLabel %s\n" % self.flowcell_name)
        hub.write("longLabel Tag sequencing, aligned %s\n" % (self.flowcell_date))
        hub.write("genomesFile genomes.txt\n")
        hub.write("email anishida@altius.org\n")
        hub.write("descriptionUrl description.html\n")
        hub.close()

    # function for creating genome.txt
    def create_genomestxt(self):
        genomefile = os.path.join(self.outdir, "genomes.txt")
        logging.info("Creating genome.txt file: %s" % genomefile)
        genomes = open(genomefile, "w")
        for hgdb in self.subtrack_sets.keys():
            genomes.write("\ngenome %s\n" % hgdb)
            genomes.write(
                "trackDb %s/trackDb.%s.%s.txt\n" % (hgdb, self.project, self.main_label)
            )
        genomes.close()

    # splits tracks up and prepares them writing
    def prepare_tracks(self):
        self.subtrack_sets = {}
        self.tracks = []

        for lane in self.data:
            logging.debug("preparing tracks for lane: " + str(lane))
            if "hgdb" not in lane:
                logging.error("Not using lane %s: no hgdb value" % lane)
                continue

            if lane["Index"] == "":
                lane["Index"] = "NoIndex"

            if lane["hgdb"] not in self.subtrack_sets:
                self.subtrack_sets[lane["hgdb"]] = []

            if lane["aligner"] == "bwa":
                track = lane.copy()
                track["strand"] = ""
                self.tracks.append(track)

            elif lane["aligner"] == "tophat":
                for strand in self.rna_strands:
                    track = lane.copy()
                    track["strand"] = strand
                    self.tracks.append(track)

        for track in self.tracks:
            hgdb = track["hgdb"]

            trackname_suffix = "L%s%s%s%sm%d" % (
                track["Lane"],
                track["Index"],
                track["SampleID"].lower(),
                track["strand"],
                self.mersize,
            )
            track["tagtrackname"] = mysql_clean(
                "%stag%s" % (self.main_label, trackname_suffix)
            )
            track["dentrackname"] = mysql_clean(
                "%sden%s" % (self.main_label, trackname_suffix)
            )

            logging.debug("tag track name: " + track["tagtrackname"])
            logging.debug("den track name: " + track["dentrackname"])

            project = track["SampleProject"]

            if self.link_dir:
                track["sampleDir"] = os.path.join(
                    "Project_%s" % project,
                    "Sample_%s" % track["SampleID"],
                    track["AlignDir"] if not poptions.pre_align_dir else "",
                )
                track["pathPrefix"] = "%s/%s" % (self.link_dir, track["sampleDir"])
            else:
                track["sampleDir"] = os.path.join(
                    self.basedir,
                    self.project_dir[project],
                    "Sample_%s" % track["SampleID"],
                    track["AlignDir"] if not poptions.pre_align_dir else "",
                )
                track["pathPrefix"] = track["sampleDir"]

            if track["aligner"] == "bwa":
                track["wigfilename"] = "%s.75_20.%s.wig" % (track["SampleName"], hgdb)
                track["bigwigfilename"] = "%s.75_20.%s.bw" % (track["SampleName"], hgdb)
                track["bamfilename"] = "%s.uniques.sorted.bam" % (track["SampleName"])
            elif track["aligner"] == "tophat":
                filename_prefix = "%s.%s.%s" % (
                    track["SampleName"],
                    track["strand"],
                    hgdb,
                )
                track["wigfilename"] = "%s.wig" % filename_prefix  # NYI
                track["bigwigfilename"] = "%s.bw" % filename_prefix
                track["bamfilename"] = "%s.bam" % filename_prefix

            # TODO: Make the RNA pipeline aware of this
            # this is to deal with the mouse with human hg19 chr11
            if hgdb == "hg19" and track["hgdb"] == "Mus_musculus":
                track["bamfilename"] = "%s_%s_L00%s.uniques.sorted.hg19.bam" % (
                    track["SampleID"],
                    track["Index"],
                    track["Lane"],
                )

            if hgdb == "hg19" and track["SampleRef"] == "Saccharomyces_cerevisiae":
                track["bamfilename"] = "%s_%s_L00%s.uniques.sorted.hg19.bam" % (
                    track["SampleID"],
                    track["Index"],
                    track["Lane"],
                )

            track["hasTags"] = False
            track["hasDensities"] = False

            if "Extra" in track and track["Extra"] is not None:
                track["Extra"] = track["Extra"].strip()

            if (
                os.path.exists(os.path.join(track["sampleDir"], track["wigfilename"]))
                and not self.bigwig
            ):
                track["hasDensities"] = True
            if (
                os.path.exists(
                    os.path.join(track["sampleDir"], track["bigwigfilename"])
                )
                and self.bigwig
            ):
                track["hasDensities"] = True
            if os.path.exists(os.path.join(track["sampleDir"], track["bamfilename"])):
                track["hasTags"] = True

            if not track["hasDensities"] or not track["hasTags"]:
                logging.error("%s does not have all files" % track["SampleID"])
                if not track["hasDensities"]:
                    logging.error("Missing densities")
                    if self.bigwig:
                        logging.error(
                            "Wanted: "
                            + os.path.join(track["sampleDir"], track["bigwigfilename"])
                        )
                    else:
                        logging.error(
                            "Wanted: "
                            + os.path.join(track["sampleDir"], track["wigfilename"])
                        )
                if not track["hasTags"]:
                    logging.error("Missing tags")
                    logging.error(
                        "Wanted: "
                        + os.path.join(track["sampleDir"], track["bamfilename"])
                    )
                logging.info("%s" % str(track))

            if track["hasDensities"] or track["hasTags"]:
                self.subtrack_sets[hgdb].append(track)

    # writes html files for individual genomes and a concatenated version with all genomes (we don't really use these)
    def create_htmls(self):
        self.html_files = {}
        masterhtmlloc = os.path.join(self.outdir, "description.html")
        masterhtml = open(masterhtmlloc, "w")
        for hgdb in self.subtrack_sets.keys():
            self.html_files[hgdb] = os.path.join(
                self.outdir, hgdb, "%s.html" % self.main_label
            )
            html = open(self.html_files[hgdb], "w")
            self.create_html(hgdb, html)
            self.create_html(hgdb, masterhtml)
            html.close()
        masterhtml.close()

    # writes the genome HTML output to file (we don't really use these)
    def create_html(self, hgdb, file):
        columns = [
            "Lane",
            "Index",
            "SampleID",
            "SampleRef",
            "CellType",
            "Assay",
            "Factors",
            "Extra",
            "wellmapping",
            "wellmapping-no-mito",
            "SPOT",
        ]
        file.write(
            "<p>Total number of lanes from this flowcell for genome %s: %d </p>\n"
            % (hgdb, len(self.subtrack_sets[hgdb]))
        )
        file.write("<table>\n")
        file.write("<thead>\n")
        file.write("<tr>\n")
        [file.write("<th>%s</th>\n" % column) for column in columns]
        file.write("</thead><tbody>\n")
        for track in self.subtrack_sets[hgdb]:
            file.write("<tr>\n")
            [file.write("<td>%s</td>\n" % track[column]) for column in columns]
            file.write("</tr>\n")
        file.write("</tbody>\n")
        file.write("</table>\n")

    def create_ras(self):
        self.ra_files = {}

        for hgdb in self.subtrack_sets.keys():
            self.create_ra(hgdb)

    # write RA / track file
    def create_ra(self, hgdb):
        logging.info("CREATING RA FOR %s" % hgdb)
        subtracks = self.subtrack_sets[hgdb]

        foldercheck(os.path.join(self.outdir, hgdb))

        self.ra_files[hgdb] = os.path.join(
            self.outdir, hgdb, "trackDb.%s.%s.txt" % (self.project, self.main_label)
        )
        ra = open(self.ra_files[hgdb], "w")

        samples = set([subtrack["SampleID"] for subtrack in subtracks])

        samples = dict()
        for subtrack in subtracks:
            if subtrack["SampleID"] not in subtrack:
                samples[subtrack["SampleID"]] = "%s %s %s %s" % (
                    subtrack["SampleID"],
                    subtrack["CellType"],
                    subtrack["Assay"],
                    subtrack["Factors"],
                )
                samples[subtrack["SampleID"]] = (
                    samples[subtrack["SampleID"]].strip().replace(" ", "_")
                )

        ra.write("track %s\n" % self.main_label)
        ra.write("compositeTrack on\n")
        ra.write("shortLabel %s\n" % self.flowcell_name)
        ra.write("longLabel Tag sequencing, aligned %s\n" % (self.flowcell_date,))
        ra.write("group %s\n" % self.project)
        ra.write("priority %s\n" % self.priority)
        ra.write("subGroup1 view Views TAG=Tags DEN=Density\n")
        ra.write(
            "subGroup2 sample Sample %s\n"
            % " ".join(
                sorted(["%s=%s" % (id, display) for id, display in samples.items()])
            )
        )
        ra.write("dimensions dimensionX=view dimensionY=sample\n")
        ra.write("sortOrder view=+ sample=+\n")
        ra.write("dragAndDrop subTracks\n")
        ra.write("type bed 3 +\n")
        ra.write("noInherit on\n\n")

        ra.write("\ttrack %stag\n" % self.main_label)
        ra.write("\tsubTrack %s\n" % self.main_label)
        ra.write("\tview TAG\n")
        ra.write("\tshortLabel Tags\n")
        ra.write("\tvisibility hide\n\n")

        for subtrack in subtracks:
            if "wellmapping-no-mito" not in subtrack:
                logging.warn(
                    "%s has no wellmapping-no-mito count" % subtrack["dentrackname"]
                )
                subtrack["wellmapping-no-mito"] = "N/A"
            if "wellmapping" not in subtrack:
                logging.warn("%s has no wellmapping count" % subtrack["dentrackname"])
                subtrack["wellmapping"] = "N/A"
            if "SPOT" not in subtrack:
                logging.warn("%s has no SPOT score" % subtrack["dentrackname"])
                subtrack["SPOT"] = "N/A"

        for subtrack in subtracks:
            ra.write("\t\ttrack %s\n" % subtrack["tagtrackname"])
            ra.write(
                "\t\tbigDataUrl %s%s/%s/%s\n"
                % (
                    self.trackhubURL,
                    self.label,
                    subtrack["sampleDir"],
                    subtrack["bamfilename"],
                )
            )
            ra.write("\t\tsubTrack %stag\n" % self.main_label)
            ra.write(
                "\t\tshortLabel %s %s:%s %s tags\n"
                % (
                    subtrack["SampleID"],
                    subtrack["Lane"],
                    subtrack["Index"],
                    subtrack["strand"],
                )
            )
            ra.write("\t\tsubGroups view=TAG sample=%s\n" % subtrack["SampleID"])
            ra.write("\t\tbamColorMode strand\n")
            ra.write(
                "\t\tlongLabel %s %s %s %s:%s %dm %s %s %s %s tags: %s (%s), spot: %s\n"
                % (
                    subtrack["CellType"],
                    subtrack["SampleID"],
                    self.flowcell_name,
                    subtrack["Lane"],
                    subtrack["Index"],
                    self.mersize,
                    subtrack["Assay"],
                    subtrack["Factors"],
                    subtrack["Extra"],
                    subtrack["strand"],
                    subtrack["wellmapping"],
                    subtrack["wellmapping-no-mito"],
                    subtrack["SPOT"],
                )
            )
            if self.paired_end:
                ra.write("\t\tpairEndsByName .\n")
            ra.write("\t\ttype bam\n\n")

        logging.info("DEN SUBTRACK GROUP")

        ra.write("\ttrack %sden\n" % self.main_label)
        ra.write("\tsubTrack %s\n" % self.main_label)
        ra.write("\tview DEN\n")
        ra.write("\tshortLabel Density\n")
        ra.write("\tvisibility full\n")
        ra.write("\tviewLimits 1:100\n")
        ra.write("\tautoScale off\n")
        ra.write("\tmaxHeightPixels 100:32:16\n\n")

        for subtrack in subtracks:
            ra.write("\t\ttrack %s\n" % subtrack["dentrackname"])
            ra.write(
                "\t\tbigDataUrl %s%s/%s/%s\n"
                % (
                    self.trackhubURL,
                    self.label,
                    subtrack["sampleDir"],
                    subtrack["bigwigfilename"],
                )
            )
            ra.write("\t\tsubTrack %sden\n" % self.main_label)
            ra.write("\t\tsubGroups view=DEN sample=%s\n" % subtrack["SampleID"])
            ra.write(
                "\t\tshortLabel %s %s:%s density\n"
                % (
                    subtrack["SampleID"],
                    subtrack["Lane"],
                    subtrack["Index"],
                )
            )
            ra.write(
                "\t\tlongLabel %s %s %s %s:%s %dm %s %s %s %s tags: %s (%s), spot: %s\n"
                % (
                    subtrack["CellType"],
                    subtrack["SampleID"],
                    self.flowcell_name,
                    subtrack["Lane"],
                    subtrack["Index"],
                    self.mersize,
                    subtrack["Assay"],
                    subtrack["Factors"],
                    subtrack["Extra"],
                    subtrack["strand"],
                    subtrack["wellmapping"],
                    subtrack["wellmapping-no-mito"],
                    subtrack["SPOT"],
                )
            )
            ra.write("\t\tgroup %s\n" % self.project)
            if self.bigwig:
                ra.write("\t\ttype bigWig\n\n")
            else:
                ra.write("\t\ttype wig 1.00 10000\n\n")

        ra.close()

    def get_sample_dir(self, lane):
        sample = lane["SampleID"]
        project = lane["SampleProject"]
        return os.path.join(self.project_dir[project], "Sample_" + sample)


class LimsQuery(object):
    def __init__(self, api_url, api_token):
        self.api_url = api_url
        self.api_token = api_token
        self.cache = dict()
        self.cache[None] = None

        self.count_types = set(["u-pf-n-mm2", "u-pf-n-mm2-mito"])

    def get(self, query):
        return self.get_by_url("%s/%s" % (self.api_url, query))

    def get_by_url(self, url):
        if url not in self.cache:
            self.cache[url] = requests.get(
                url, headers={"Authorization": "Token %s" % self.api_token}
            ).json()
        return self.cache[url]

    def get_all(self, query):
        data = self.get(query)
        results = data["results"]
        while data["next"] is not None:
            data = self.get_by_url(data["next"])
            results += data["results"]
        return results

    def get_counttype_by_codename(self, codename):
        return self.get_all("flowcell_lane_count_type/?codename=%s" % codename)[0]

    def get_counts_for_alignment(self, alignment):
        counts = dict()
        for type in self.count_types:
            type_id = self.get_counttype_by_codename(type)["id"]
            count_vals = self.get_all(
                "flowcell_lane_count/?alignment=%s&count_type=%d" % (alignment, type_id)
            )
            if count_vals:
                counts[type] = count_vals[0]["count"]

        # Check to see if we got all the types we wanted
        for count in self.count_types:
            if count not in counts:
                logging.warn(
                    "Could not fetch count %s for alignment: %s" % (count, alignment)
                )

        return counts

    def get_rna_metrics_for_alignment(self, alignment):
        results = self.get("rna_alignment_metrics/?alignment=%s" % alignment)
        if not results["results"]:
            logging.warn("Could not fetch RNA metrics for alignment: %s" % alignment)
            return None
        return results["results"][0]

    def get_alignment(self, id):
        return self.get("flowcell_lane_alignment/%s/" % id)

    def get_spot_for_alignment(self, alignment):
        # TODO: This assumes one spot per alignment.
        results = self.get("flowcell_lane_spot/?alignment=%s" % alignment)
        if not results["results"]:
            return None
        return results["results"][0]


def get_alignment_data(library, alignment, lims):
    # This is mainly a shim.

    logging.debug("Fetching data for library: %s" % library)

    d = dict()
    d["project"] = library["project"]
    d["hgdb"] = alignment["genome_index"]
    d["aligner"] = alignment["aligner"]
    d["SampleName"] = alignment["sample_name"]
    d["AlignDir"] = alignment["align_dir"]
    d["Index"] = library["barcode_index"]
    d["SampleID"] = library["samplesheet_name"]
    # cell_type included for backwards compatibility with older processing files (before Feb 2016)
    d["CellType"] = library.get("sample_taxonomy") or library.get("cell_type")
    d["Assay"] = library["assay"]
    d["Lane"] = library["lane"]
    d["SampleProject"] = library["project"]

    lims_lane = lims.get("flowcell_lane/%s" % library["id"])
    lims_sample = lims.get_by_url(lims_lane["sample"])

    d["failed_lane"] = lims_lane["failed"]
    if d["failed_lane"]:
        logging.warn("Lane marked as failed, not using: %s" % library["id"])
        return d

    if d["aligner"] == "bwa":
        lims_counts = lims.get_counts_for_alignment(alignment["id"])
        d["wellmapping"] = lims_counts.get("u-pf-n-mm2", None)
        d["wellmapping-no-mito"] = lims_counts.get("u-pf-n-mm2-mito", None)

    # RNA doesn't have u-pf-no-mito counts, so we set those properties from the rna metrics
    elif d["aligner"] == "tophat":
        r = lims.get_rna_metrics_for_alignment(alignment["id"])
        if r is not None:
            # Subtract off ribosomal RNA
            d["wellmapping"] = int(r["mapped_reads"])
            d["wellmapping-no-mito"] = int(
                int(r["mapped_reads"]) * (1 - (float(r["percent_chrM"]) / 100.0))
            )

    d["Extra"] = lims_lane["extra"]
    d["SampleRef"] = ""  # NYI
    if lims_sample is not None:
        d["Factors"] = ", ".join(
            [
                lims.get_by_url(factor)["display_name"]
                for factor in lims_sample["factors"]
            ]
        )
    else:
        d["Factors"] = None

    lims_spot = lims.get_spot_for_alignment(alignment["id"])
    d["SPOT"] = lims_spot["spot_score"] if lims_spot else "N/A"

    return d


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

    data = json.loads(open(poptions.process_config, "r").read())
    trackhubconfig = poptions.trackhub_config

    # projects = [d["code_name"] for d in data["projects"]]

    # get basedir
    basedir = data["alignment_group"]["directory"]

    # for old flowcells where the LIMS hasn't been updated, update old file prefix to new prefix
    # (this might be fixed now?)
    if re.search("monarch", basedir):
        fc_loc = "/net/seq/data/flowcells/"
        fc_dirname = basedir.replace("/net/monarch/vol2/tag/stamlab/flowcells/", "")
        basedir = fc_loc + fc_dirname

    # fetch paired endedness
    paired_end = data["flowcell"]["paired_end"]
    label = data["alignment_group"]["label"]
    mersize = data["flowcell"]["read_length"]
    date = data["alignment_group"]["label"].split("_")[1]

    # get browsersheet information
    lims = LimsQuery(poptions.api_url, poptions.api_token)

    load_groups = dict()

    # find projects
    for lib in data["libraries"]:
        for a in lib["alignments"]:
            align_data = get_alignment_data(lib, a, lims)
            if not align_data["failed_lane"]:
                p = align_data["project"]
                if p not in load_groups:
                    load_groups[p] = []
                load_groups[p].append(align_data)

    for project in load_groups.keys():
        lane_group = load_groups[project]
        logging.info("the basedirectory is: %s" % basedir)
        outdir = os.path.join(basedir, "browser-load-%s" % project)
        loader = MakeBrowserload(
            lane_group,
            trackhubconfig,
            basedir,
            outdir,
            mersize,
            poptions.priority,
            paired_end,
            project,
            label,
            date,
        )
        loader.load()


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
