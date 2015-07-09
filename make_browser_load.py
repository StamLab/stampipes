#!/usr/bin/env python
import json
import os
import sys
import argparse
import logging
import re
import copy
import requests

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
logging.getLogger("requests").setLevel(logging.WARNING)

options = {
    "quiet": False,
    "debug": False,
    "process_config": "processing.json",
    "priority": None,
    "api_url": os.getenv('LIMS_API_URL'),
    "api_token": os.getenv('LIMS_API_TOKEN'),
}

util_log = logging.getLogger("StamPy.util")
def foldercheck(*args):
    """Checks to see if the folders exist, creates them if they are not."""

    for folder in args:
        if not os.path.isdir(folder):
            try:
                os.mkdir(folder)
                util_log.info("Created folder: %s" % folder)
            except OSError, x:
                util_log.error("ERROR: Could not create directory: %s" % folder)
                util_log.warn("Please make sure all nonexistant parent directories have been created.")
                sys.exit(0)

def mysql_clean(input):
    # Mysql names can contain only 0-9, a-z, A-Z, _, or $
    # So we replace all other characters with an underscore,
    # after removing leading/trailing whitespace.
    output = re.sub("[^\w$]", "_", input.strip())
    return output

def parser_setup():
    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-j", "--json-config", dest="process_config",
        help="The process config to work off of.")
    parser.add_argument("-p", "--priority", dest="priority", required=True,
        help="The priority of this flowcell")
    parser.add_argument("-n", "--post-aggregation", dest="post_aggregation", action="store_true",
        help="Pass this option to indicate that the directory structure for this flowcell is post-aggregation")

    parser.set_defaults( **options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

class MakeBrowserload(object):
    genome_organisms = {
      "hg19": "human",
      "rn5": "rat",
      "mm9": "mouse",
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

    rna_strands = [ "all", "pos", "neg" ]

    #def __init__(self, browserconfig, browsersheet, basedir, outdir, priority, paired_end, project, project_dir = "",
        #maintrackname = None, bigwig = True, date = None):
    def __init__(self, group_data, browserconfig, basedir, outdir, priority, paired_end, project, date):

        self.basedir = basedir
        self.flowcell_date = date
        self.outdir = outdir
        self.mersize = 36
        self.win=75
        self.binI=20
        self.priority = priority
        self.paired_end = paired_end
        self.bigwig = True
        self.projects = project.split(",")
        self.date = date
        self.project_dirs = {}
        self.data = group_data
        project_dir = ""

        if len(self.projects) == 1 and project_dir:
            self.project_dirs[project] = project_dir
            logging.info("Using project dir: %s" % self.project_dirs[project])
        else:
            for project in self.projects:
                self.project_dirs[project] = os.path.join(self.basedir, "Project_" + project)

        self.load_config(browserconfig)

    def load_config(self, browserconfig):
       import ConfigParser
       Config = ConfigParser.ConfigParser()
       Config.read(browserconfig)
       self.server = Config.get("browser", "server")
       self.browser_url = Config.get("browser", "browser_url")
       self.flowcell_link_folder = Config.get("browser", "flowcell_link_folder")
       self.track_basedir = Config.get("browser", "track_basedir")
       self.browser_excludes_file = Config.get("browser", "browser_excludes_file")
       self.group = Config.get('browser', 'browser_group')
       self.file_label = Config.get('browser', 'file_label')

    def load(self):
        #self.browsersheet = SampleSheet(file=self.browsersheet_file)
        self.basedir_name = os.path.basename(self.basedir)
        foldercheck(self.outdir)

        #if self.maintrackname:
        if False:
           self.main_label = "%s%son%s" % (self.file_label, self.maintrackname, self.date)
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

            self.main_label = "%s%son%s" % (self.file_label, self.flowcell_name, self.flowcell_date)

        logging.info("Main track name: %s" % self.main_label)

        self.excludes_file = os.path.join(self.outdir, "excludes.%s" % self.main_label)

        if self.flowcell_link_folder:
            logging.debug("link folder: " + self.flowcell_link_folder + " base folder: " + self.basedir_name)
            self.link_dir = os.path.join(self.flowcell_link_folder, self.basedir_name)
        else:
            self.link_dir = ""

        self.prepare_tracks()
        logging.info("Main label: %s" % self.main_label)

        self.create_ras()
        self.create_htmls()
        self.create_commands()
        self.create_excludes()

    def prepare_tracks(self):
        """Splits the tracks up and makes changes to the data to make it easier for later on."""

        self.subtrack_sets = {}
        self.tracks = []

        for lane in self.data:

            logging.debug("preparing tracks for lane: " + str(lane))
            if not "hgdb" in lane:
                logging.error("Not using lane %s: no hgdb value" % lane )
                continue

            if lane["Index"] == "":
                lane["Index"] = "NoIndex"

            if not lane["hgdb"] in self.subtrack_sets:
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

            trackname_suffix = "L%s%s%s%sm%d" % (track["Lane"], track["Index"], track["SampleID"].lower(), track["strand"], self.mersize)
            track["tagtrackname"] = mysql_clean("%stag%s" % (self.main_label, trackname_suffix))
            track["dentrackname"] = mysql_clean("%sden%s" % (self.main_label, trackname_suffix))

            logging.debug("tag track name: " + track["tagtrackname"])
            logging.debug("den track name: " + track["dentrackname"])

            project = track["SampleProject"]

            if self.link_dir:
                track["sampleDir"] = os.path.join("Project_%s" % project, "Sample_%s" % track["SampleID"])
                track["pathPrefix"] = "%s/%s" % (self.link_dir, track["sampleDir"])
            else:
                track["sampleDir"] = os.path.join(self.basedir, self.project_dir[project], "Sample_%s" % track["SampleID"])
                track["pathPrefix"] = track["sampleDir"]
            if poptions.post_aggregation: 
                track["sampleDir"] = os.path.join(track["sampleDir"], track['AlignDir'])

            if lane["aligner"] == "bwa":
                track["wigfilename"]    = "%s.75_20.%s.wig"       % (track["SampleName"], hgdb)
                track["bigwigfilename"] = "%s.75_20.%s.bw"        % (track["SampleName"], hgdb)
                track["bamfilename"]    = "%s.uniques.sorted.bam" % (track["SampleName"])
            elif lane["aligner"] == "tophat":
                filename_prefix = "%s.%s.%s"     % (track["SampleName"], track["strand"], hgdb)
                track["wigfilename"]    = "%s.wig" % filename_prefix  # NYI
                track["bigwigfilename"] = "%s.bw"  % filename_prefix
                track["bamfilename"]    = "%s.bam" % filename_prefix


            # TODO: Make the RNA pipeline aware of this
            # this is to deal with the mouse with human hg19 chr11
            if( hgdb == "hg19" and track["hgdb"] == "Mus_musculus" ):
                track["bamfilename"] = "%s_%s_L00%s.uniques.sorted.hg19.bam" % (track["SampleID"], track["Index"], track["Lane"])

            if( hgdb == "hg19" and track["SampleRef"] == "Saccharomyces_cerevisiae" ):
                track["bamfilename"] = "%s_%s_L00%s.uniques.sorted.hg19.bam" % (track["SampleID"], track["Index"], track["Lane"])

            track["hasTags"] = False
            track["hasDensities"] = False

            track["Extra"] = track["Extra"].strip()

            if os.path.exists(os.path.join(track["sampleDir"], track["wigfilename"])) and not self.bigwig:
                track["hasDensities"] = True
            if os.path.exists(os.path.join(track["sampleDir"], track["bigwigfilename"])) and self.bigwig:
                track["hasDensities"] = True
            if os.path.exists(os.path.join(track["sampleDir"], track["bamfilename"])):
                track["hasTags"] = True

            if not track["hasDensities"] or not track["hasTags"]:
                logging.error("%s does not have all files" % track["SampleID"])
                if not track["hasDensities"]:
                    logging.error( "Missing densities" )
                    if self.bigwig:
                        logging.error("Wanted: " + os.path.join(track["sampleDir"], track["bigwigfilename"]))
                    else:
                        logging.error("Wanted: " + os.path.join(track["sampleDir"], track["wigfilename"]))
                if not track["hasTags"]:
                    logging.error("Missing tags")
                    logging.error("Wanted: " + os.path.join(track["sampleDir"], track["bamfilename"]))
                logging.info("%s" % str(track))

            if track["hasDensities"] or track["hasTags"]:
                self.subtrack_sets[hgdb].append(track)

    def create_htmls(self):
        self.html_files = {}

        for hgdb, subtracks in self.subtrack_sets.items():
            self.create_html(hgdb)

    def create_html(self, hgdb):
        self.html_files[hgdb] = os.path.join(self.outdir, hgdb, "%s.html" % self.main_label)

        html = open( self.html_files[hgdb], 'w')

        columns = ["Lane", "Index", "SampleID", "SampleRef", "CellType", "Assay", "Factors", "Extra",
            "wellmapping", "wellmapping-no-mito", "SPOT"]

        html.write("<p>Total number of lanes from this flowcell for this genome: %d </p>\n" % len(self.subtrack_sets[hgdb]))

        html.write("<table>\n")
        html.write("<thead>\n")
        html.write("<tr>\n")
        [html.write("<th>%s</th>\n" % column) for column in columns]
        html.write("</thead><tbody>\n")

        for track in self.subtrack_sets[hgdb]:
            html.write("<tr>\n")
            [html.write("<td>%s</td>\n" % track[column]) for column in columns]
            html.write("</tr>\n")

        html.write("</tbody>\n")
        html.write("</table>\n")

        html.close()

    def create_ras(self):
        self.ra_files = {}

        for hgdb, subtracks in self.subtrack_sets.items():
            self.create_ra(hgdb)

    def create_commands(self):
        makefile = os.path.join(self.outdir, "make.%s.doc" % self.main_label)
        logging.info("Makefile: %s" % makefile)
        commands = open( makefile, 'w')

        commands.write("# %s\n" % makefile)
        commands.write("# %s\n\n" % ", ".join(self.subtrack_sets.keys()))

        if self.link_dir:
            commands.write("ln -s %s %s\n\n" % (self.basedir, self.link_dir))

        for hgdb, subtracks in self.subtrack_sets.items():
            self.create_genome_commands( hgdb, commands)

        commands.write("\ncat %s >> %s\n" % (self.excludes_file, self.browser_excludes_file))

        commands.close()

    def create_subtrack_commands(self, subtrack, commandsout):
        if subtrack["hasDensities"] and not self.bigwig:
            commandsout.write("hgLoadWiggle -pathPrefix=%s %s %s %s/%s\n" % (
            subtrack["pathPrefix"], subtrack["hgdb"], subtrack["dentrackname"], subtrack["pathPrefix"], subtrack["wigfilename"]))
        # hgLoadWiggle -pathPrefix=/usr/local/UW/flowcell-density/FCB0BLA_110620_tag/005 hg19 STAM_FCB0BLA_110620_IT_DEN_L005_6_DS18466_36_DNaseI /usr/local/UW/flowcell-density/FCB0BLA_110620_tag/005/FCB0BLA_lane6_75_20.wig

        if subtrack["hasDensities"] and self.bigwig:
            commandsout.write("hgBbiDbLink %s %s %s/%s\n" % (subtrack["hgdb"], subtrack["dentrackname"], subtrack["pathPrefix"], subtrack["bigwigfilename"]))

        if subtrack["hasTags"]:
            hgsqlcommand = "hgsql %s -e '" % subtrack["hgdb"]
            hgsqlcommand += "drop table if exists %s; " % subtrack["tagtrackname"]
            hgsqlcommand += "create table %s (filename varchar(255) not null); " % subtrack["tagtrackname"]
            hgsqlcommand += "insert into %s values " % subtrack["tagtrackname"]
            hgsqlcommand += "(\"%s/%s\");'\n" % (subtrack["pathPrefix"], subtrack["bamfilename"])

            commandsout.write(hgsqlcommand)

#ln -s $datafile bam-links/Rudensky/Rudensky_bams/$data.bam
#ln -s $indexfile bam-links/Rudensky/Rudensky_bams/$data.bam.bai
#hgsql $forg -e 'drop table if exists $trackType; create table
#$trackType (fileName varchar(255) not null); insert into $trackType
#values (\"/usr/local/UW/bam-links/Rudensky/Rudensky_bams/$data.bam\");'"

    def create_genome_commands(self, hgdb, commandsout):
        if not hgdb in self.genome_organisms:
            logging.error(hgdb + " not in " + str(self.genome_organisms))
            commandsout.write("\n ERROR: no " + hgdb + " genome\n")
            return

        organism = self.genome_organisms[hgdb]

        commandsout.write("# Creating commands for %s\n\n" % hgdb)
        for subtrack in self.subtrack_sets[hgdb]:
            self.create_subtrack_commands(subtrack, commandsout)

        commandsout.write("\ncp %s %s/%s/%s\n" % (self.html_files[hgdb], self.track_basedir, organism, hgdb))
        commandsout.write("cp %s %s/%s/%s\n\n" % (self.ra_files[hgdb], self.track_basedir, organism, hgdb))
        commandsout.write('# add line "include trackDb.%s.%s.ra" to %s/%s/%s/trackDb.%s.ra\n\n' % (self.file_label, self.main_label, self.track_basedir, organism, hgdb, self.file_label))

    def create_excludes(self):
        excludes = open( self.excludes_file, 'w')

        for subtrack in self.tracks:
            for suffix in ["frm", "MYD", "MYI"]:
                logging.debug( "subtrack contents: " + str(subtrack))
                excludes.write("%s.%s\n" % (subtrack["tagtrackname"], suffix))
                excludes.write("%s.%s\n" % (subtrack["dentrackname"], suffix))

        excludes.close()

    def create_ra(self, hgdb):
        logging.info("CREATING RA FOR %s" % hgdb)
        subtracks = self.subtrack_sets[hgdb]

        foldercheck(os.path.join(self.outdir, hgdb))

        self.ra_files[hgdb] = os.path.join(self.outdir, hgdb, "trackDb.%s.%s.ra" % (self.file_label, self.main_label))
        ra = open( self.ra_files[hgdb], 'w' )

        samples = set([subtrack["SampleID"] for subtrack in subtracks])

        samples = dict()
        for subtrack in subtracks:
            if not subtrack["SampleID"] in subtrack:
                samples[subtrack["SampleID"]] = "%s %s %s %s" % (subtrack["SampleID"], subtrack["CellType"], subtrack["Assay"], subtrack["Factors"])
                samples[subtrack["SampleID"]] = samples[subtrack["SampleID"]].strip().replace(" ", "_")

        ra.write("track %s\n" % self.main_label)
        ra.write("compositeTrack on\n")
        ra.write("shortLabel %s\n" % self.flowcell_name)
        ra.write("longLabel Tag sequencing, aligned %s\n" % (self.flowcell_date,))
        ra.write("group %s\n" % self.group)
        ra.write("priority %s\n" % self.priority)
        ra.write("subGroup1 view Views TAG=Tags DEN=Density\n")
        ra.write("subGroup2 sample Sample %s\n" % " ".join(sorted(['%s=%s' % (id, display) for id, display in samples.iteritems()])))
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
            if not "wellmapping-no-mito" in subtrack:
                logging.warn("%s has no wellmapping-no-mito count" % subtrack["dentrackname"] )
                subtrack["wellmapping-no-mito"] = "N/A"
            if not "wellmapping" in subtrack:
                logging.warn("%s has no wellmapping count" % subtrack["dentrackname"] )
                subtrack["wellmapping"] = "N/A"
            if not "SPOT" in subtrack:
                logging.warn("%s has no SPOT score" % subtrack["dentrackname"] )
                subtrack["SPOT"] = "N/A";

        #track STAM_FC630D3_110711_IT_TAG_L5_DS18900_36_
        #       subTrack STAM_FC630D3_110711_IT_TAG
        #       subGroups view=TAG
        #       shortLabel DS18900 5 tags
        #       longLabel D_GM12864_DS18900_I_FC630D3_5_36m | GM12864  [Expansion] DS18900 110711 FC630D3 L5 - n tags = 19112570 (19108346) : ptih: 0.5589
        #       group illumina-raw
        #       type bed 6

        for subtrack in subtracks:
            ra.write("\t\ttrack %s\n" % subtrack["tagtrackname"])
            ra.write("\t\tsubTrack %stag\n" % self.main_label)
            ra.write("\t\tshortLabel %s %s:%s %s tags\n" % (subtrack["SampleID"], subtrack["Lane"], subtrack["Index"], subtrack["strand"]))
            ra.write("\t\tsubGroups view=TAG sample=%s\n" % subtrack["SampleID"])
            ra.write("\t\tbamColorMode strand\n")
            ra.write("\t\tlongLabel %s %s %s %s:%s %dm %s %s %s %s tags: %s (%s), spot: %s\n" % (
                subtrack["CellType"], subtrack["SampleID"], self.flowcell_name, subtrack["Lane"],
                subtrack["Index"], self.mersize, subtrack["Assay"], subtrack["Factors"], subtrack["Extra"], subtrack["strand"], subtrack["wellmapping"],
                subtrack["wellmapping-no-mito"], subtrack["SPOT"]))
            if self.paired_end:
                ra.write("\t\tpairEndsByName .\n")
            ra.write("\t\ttype bam\n\n")

        logging.info("DEN SUBTRACK GROUP")

        ra.write( "\ttrack %sden\n" % self.main_label)
        ra.write( "\tsubTrack %s\n" % self.main_label)
#        track STAM_FC630D3_110711_IT_DEN
 #       subTrack STAM_FC630D3_110711_IT

        ra.write("\tview DEN\n")
        ra.write("\tshortLabel Density\n")
        ra.write("\tvisibility full\n")
        ra.write("\tviewLimits 1:100\n")
        ra.write("\tautoScale off\n")
        ra.write("\tmaxHeightPixels 100:32:16\n\n")

        #  track STAM_FC62J4G_101107_IT_DEN_L5_DS13475_36_DNaseI
        #        subTrack STAM_FC62J4G_101107_IT_DEN
        #        subGroups view=DEN
        #        shortLabel DS13475 5 density
        #        longLabel D_HUVEC_DS13475_I_FC62J4G_5_36m | HUVEC DNaseI [Expansion] DS13475 101107 FC62J4G L5 - n tags = 22137427 (20360302) : ptih: 0.2946 density
        #        group illumina-raw
        #        type wig 1.00 10000

        for subtrack in subtracks:
            ra.write("\t\ttrack %s\n" % subtrack["dentrackname"])
            ra.write("\t\tsubTrack %sden\n" % self.main_label)
            ra.write("\t\tsubGroups view=DEN sample=%s\n" % subtrack["SampleID"])
            ra.write("\t\tshortLabel %s %s:%s density\n" % (subtrack["SampleID"], subtrack["Lane"], subtrack["Index"],))
            ra.write("\t\tlongLabel %s %s %s %s:%s %dm %s %s %s %s tags: %s (%s), spot: %s\n" % (
                subtrack["CellType"], subtrack["SampleID"], self.flowcell_name, subtrack["Lane"],
                subtrack["Index"], self.mersize, subtrack["Assay"], subtrack["Factors"], subtrack["Extra"], subtrack["strand"], subtrack["wellmapping"],
                subtrack["wellmapping-no-mito"], subtrack["SPOT"]))
            ra.write("\t\tgroup %s\n" % self.group)
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

        self.count_types = set(['u-pf-n-mm2', 'u-pf-n-mm2-mito'])

    def get(self, query):
        return self.get_by_url( "%s/%s" % ( self.api_url, query ) )

    def get_by_url(self, url):
        if not url in self.cache:
            #print url
            self.cache[url] = requests.get(url, headers={'Authorization': "Token %s" % self.api_token}).json()
        return self.cache[url]

    def get_counts_for_alignment(self, alignment):
        counts = dict()
        page = 1
        while len(counts.keys()) < len(self.count_types):
            result = self.get("flowcell_lane_count/?alignment=%s&page=%s" % (alignment, page ))
            if not 'results' in result:
                break
            for count in result['results']:
                count_name = count['count_type_name']
                if count_name in self.count_types:
                    counts[count_name] = count['count']
            page += 1

        # Check to see if we got all the types we wanted
        for count  in self.count_types:
            if count not in counts:
                logging.warn("Could not fetch count %s for alignment: %s" % (count, alignment))

        return counts

    def get_rna_metrics_for_alignment(self, alignment):
        results = self.get("rna_alignment_metrics/?alignment=%s" % alignment)
        if not results['results']:
            logging.warn("Could not fetch RNA metrics for alignment: %s" % alignment)
            return None
        return results['results'][0]

    def get_alignment(self, id):
        return self.get("flowcell_lane_alignment/%s/" % id)

    def get_spot_for_alignment(self, alignment):
        # TODO: This assumes one spot per alignment.
        results = self.get("flowcell_lane_spot/?alignment=%s" % alignment)
        if not results['results']:
            return None
        return results['results'][0]


def get_alignment_data(library, alignment, lims):
    # This is mainly a shim.

    logging.debug("Fetching data for library: %s" % library)

    d = dict()
    d['project']       = library['project']
    d['hgdb']          = alignment['genome_index']
    d['aligner']       = alignment['aligner']
    d['SampleName']    = alignment['sample_name']
    d['AlignDir']      = alignment['align_dir']
    d['Index']         = library['barcode_index']
    d['SampleID']      = library['samplesheet_name']
    d['CellType']      = library['cell_type']
    d['Assay']         = library['assay']
    d['Lane']          = library['lane']
    d['SampleProject'] = library['project']

    lims_lane = lims.get("flowcell_lane/%s" % library['id'])
    lims_sample = lims.get_by_url( lims_lane['sample'] )
    lims_library = lims.get_by_url( lims_lane['library'] )

    d['failed_lane'] = lims_lane['failed']
    if d['failed_lane']:
        logging.warn("Lane marked as failed, not using: %s" % library['id'])
        return d

    if d['aligner'] == 'bwa':
        lims_counts = lims.get_counts_for_alignment(alignment['id'])
        d['wellmapping']         = lims_counts.get('u-pf-n-mm2', None)
        d['wellmapping-no-mito'] = lims_counts.get('u-pf-n-mm2-mito', None)

    # RNA doesn't have u-pf-no-mito counts
    # So we set those properties from the rna metrics
    elif d['aligner'] == 'tophat':
        r = lims.get_rna_metrics_for_alignment(alignment['id'])
        if r is not None:
            # Subtract off ribosomal RNA
            d['wellmapping']         = int(r['mapped_reads'])
            d['wellmapping-no-mito'] = int(int(r['mapped_reads']) * (1 - (float(r['percent_chrM']) / 100.0)))

    d['Extra']         = lims_lane['extra']
    d['SampleRef']     = ""  #NYI
    if lims_sample is not None:
        d['Factors']       = ", ".join([ lims.get_by_url(factor)['display_name'] for factor in lims_sample['factors'] ])
    else:
        d['Factors'] = None


    lims_spot = lims.get_spot_for_alignment(alignment['id'])
    d['SPOT'] = lims_spot['spot_score'] if lims_spot else "N/A"

    return d

def main(args = sys.argv):
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

    data = json.loads(open(poptions.process_config, 'r').read())

    projects = [ d['code_name'] for d in data['projects'] ]

    # get basedir
    basedir = data['alignment_group']['directory']
    # Fetch paired endedness?
    paired_end = data['flowcell']['paired_end']
    date = data['alignment_group']['label'].split('_')[1]

    # get browsersheet information!

    lims = LimsQuery( poptions.api_url, poptions.api_token )

    browsers = set()
    load_groups = dict()
    for l in data['libraries']:
        for a in l['alignments']:
            if a['browsers']:    # Only process alignments that map to a browser
                align_data = get_alignment_data(l, a, lims)
                if not align_data['failed_lane']:
                    for b in a['browsers']:
                        browsers.add(b)
                        key = ( align_data['project'], b )
                        if not key in load_groups:
                            load_groups[key] = []
                        load_groups[key].append( align_data )


    for group_key in load_groups.keys():
        (project, browser) = group_key

        lane_group = load_groups[group_key]

        browserconfig = os.path.join( os.getenv("STAMPIPES"), "config", "ucsc_browser", "%s-%s.config" % (browser, project) )
        print browserconfig

        outdir = os.path.join( basedir, "browser-load-%s-%s" % (project, browser))

        loader = MakeBrowserload(lane_group, browserconfig, basedir, outdir, poptions.priority, paired_end, project, date)
        loader.load()


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
