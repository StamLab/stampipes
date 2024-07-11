#!/usr/bin/env python

"""This script is meant to be run in the background.  It will keep checking for
new sequencer folders and for sequencer folders that have finished copying.

TODO: Generalize some of these functions into a library.
TODO: This might be more robust as a cron job."""

import json
import logging
import os
import re
import smtplib
import time
import xml
from email.mime.text import MIMEText
from xml.dom import minidom

"""
{
  "base_folders": {
    "server": "/path/to/sequencing/drop"
  },
  "notification_emails": {
    "Name": "sendto@email.com",
  },
  "send_from": "sendfrom@email.com"
}
"""
config_filename = os.path.join(
    os.environ["STAMPIPES_DATA"], "flowcell_notify_config.json"
)

# This file will be checked for in the folders
# will need to append it with the read #
copy_complete_filename = "Basecalling_Netcopy_complete_Read%d.txt"

# The name of the RunInfo.xml file in the base of the flowcell directory
runinfo_filename = "RunInfo.xml"

# Format of folders: 090810_SOLEXA-1GA-1_0016_FC82IU
folder_pattern_ga = re.compile(
    r"(?P<date>\d{6})_SOLEXA-1GA-[12]_\d{4,5}_FC(?P<flowcell>[A-Z0-9]{5})"
)
# 140703_SN373_0524_BC6TATACXX
# 140710_D00453_0080_AC5PBPANXX
folder_pattern_hiseq = re.compile(
    r"(?P<date>\d{6})_(?P<sequencer_id>(SN|D)\d+)_[0-9]+_(A|B)(?P<flowcell>[A-Z0-9]{5})[A-Z]{2}XX"
)
# 140808_NS500372_0009_AH115HBGXX
folder_pattern_nextseq = re.compile(
    r"(?P<date>\d{6})_(?P<sequencer_id>NS500\d+)_[0-9]+_A(?P<flowcell>[A-Z0-9]{5})[A-Z]{2}XX"
)

# To use with datetime.strptime() to translate folder dates
folder_datepattern = "%y%m%d"

# Waiting period, in seconds
wait = 60 * 60  # every hour

config_json = open(config_filename)

config = json.loads(config_json.read())

config_json.close()

# Base folders to check for sequencing run folders
base_folders = config["base_folders"]

# WARNING: CURRENTLY ONLY STAMLAB.ORG EMAILS WORK

notification_emails = config["notification_emails"].values()

# Contents of new folder emails
notif_new_title = "New sequencer folder: %(flowcell)s (%(server)s)"
notif_new_body = "A new folder has started copying:\n\n%(folder)s"
# Contents of finished copy emails
notif_copy_title = "Finished sequencing copy: %(flowcell)s (%(server)s)"
notif_copy_body = "A sequencer run has finished copying:\n\n%(folder)s"

# The folders to be checked will be kept track of here
check_folders = dict()
flowcell_reads = dict()


def get_sequencer_folders():
    """Gets all the sequence folders in the base folders."""
    folder_list = list()

    for base_folder in base_folders.values():
        [
            folder_list.append("%s/%s" % (base_folder, sequencer_folder))
            for sequencer_folder in os.listdir(base_folder)
            if os.path.isdir("%s/%s" % (base_folder, sequencer_folder))
            and (
                folder_pattern_ga.match(sequencer_folder)
                or folder_pattern_hiseq.match(sequencer_folder)
                or folder_pattern_nextseq.match(sequencer_folder)
            )
        ]
    return folder_list


def get_folder_reads(sequencer_folder):
    runinfo_file = os.path.join(sequencer_folder, runinfo_filename)
    try:
        runinfodoc = minidom.parse(runinfo_file)
        return len(runinfodoc.getElementsByTagName("Read"))
    except IOError:
        logging.info("Could not read %s" % runinfo_file)
        return None
    except xml.parsers.expat.ExpatError:
        logging.info("%s is malformatted" % runinfo_file)
        return None


def check_copy(sequencer_folder):
    """Checks to see if the given copy filename is present in the sequencer folder"""

    if flowcell_reads[sequencer_folder]:
        # copy_filename = copy_complete_filename % flowcell_reads[sequencer_folder]
        return os.path.exists("%s/%s" % (sequencer_folder, "CopyComplete.txt"))
    else:
        return False


def load_folders():
    """This function loads the initial folder states."""
    logging.info("Loading folders")
    for sequencer_folder in get_sequencer_folders():
        flowcell_reads[sequencer_folder] = get_folder_reads(sequencer_folder)

        check_folders[sequencer_folder] = check_copy(sequencer_folder)

        if flowcell_reads[sequencer_folder]:
            logging.info(
                "Initial state of %s: %s"
                % (sequencer_folder, str(check_folders[sequencer_folder]))
            )
        else:
            logging.info("Initial state of %s: does not have reads" % sequencer_folder)


def check_folder(sequencer_folder):
    """Check a sequencer folder and notify for changes"""
    logging.debug("Checking folder: %s" % sequencer_folder)

    if sequencer_folder in check_folders:
        if not check_folders[sequencer_folder] and check_copy(sequencer_folder):
            logging.info("Folder finished copying: %s" % sequencer_folder)
            notify_copy(sequencer_folder)
            check_folders[sequencer_folder] = True
    else:
        logging.info("New folder: %s" % sequencer_folder)
        check_folders[sequencer_folder] = False
        flowcell_reads[sequencer_folder] = get_folder_reads(sequencer_folder)
        logging.debug("Number of reads: %s" % str(flowcell_reads[sequencer_folder]))
        notify_new(sequencer_folder)


def run_check():
    """Checks for newly created folders and newly copied folders.  Notifies if they are found."""
    try:
        folders = get_sequencer_folders()
    except OSError:
        logging.error("Error getting sequencer folder; skipping this cycle.")
        return None

    logging.debug("Running check")

    # delete folders that don't exist anymore from checking
    for sequencer_folder in folders:
        if not os.path.exists(sequencer_folder):
            logging.info("Deleting folder: %s" % sequencer_folder)
            del check_folders[sequencer_folder]

    # check each folder for being new or being copied
    for sequencer_folder in folders:
        check_folder(sequencer_folder)


def get_folder_info(sequencer_folder):
    """Given a sequencer folder name, tries to find the server and flowcell name"""
    info = dict()

    for server in base_folders.keys():
        if sequencer_folder.find(server) > 0:
            info["server"] = server

    if "server" not in info:
        info["server"] = "UNKNOWN"

    match = folder_pattern_ga.search(sequencer_folder)

    if not match:
        match = folder_pattern_hiseq.search(sequencer_folder)

    if not match:
        match = folder_pattern_nextseq.search(sequencer_folder)

    info["flowcell"] = "FC" + match.groupdict()["flowcell"]
    info["folder"] = sequencer_folder
    return info


def notify_new(sequencer_folder):
    """Notify the list of emails of a finished copying folder."""
    info = get_folder_info(sequencer_folder)
    send_email(notif_new_title % info, notif_new_body % info, notification_emails)


def notify_copy(sequencer_folder):
    """Notify the list of emails of a newly created folder."""
    info = get_folder_info(sequencer_folder)
    send_email(notif_copy_title % info, notif_copy_body % info, notification_emails)


def send_email(subject, body, emails):
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = "illumina@stamlab.org"
    msg["To"] = ", ".join(emails)

    s = smtplib.SMTP("localhost")
    s.sendmail("illumina@stamlab.org", emails, msg.as_string())
    s.quit()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # loads what's presently there so we don't notify about it
    load_folders()

    while True:
        # because somebody running the script already has checked the state of things beforehand,
        # wait period goes first
        logging.info("Waiting %d seconds before next check" % wait)
        time.sleep(wait)
        run_check()
