"""
For a tab delineated file with rows of:

FILENAME\tMD5SUM

Check that each file matches the needed MD5SUM
"""

import logging
import subprocess
import sys

infile = sys.argv[1]


def check_md5sum(filename, md5sum):
    logging.debug("Checking file %s matches %s" % (filename, md5sum))
    current_md5sum = subprocess.check_output(
        ["md5sum", filename], stderr=subprocess.STDOUT, universal_newlines=True
    ).split()[0]

    match = md5sum == current_md5sum

    if not match:
        logging.error(
            "md5sum for file %s does not match: %s recorded, %s as exists"
            % (filename, md5sum, current_md5sum)
        )

    return match


MATCHING = True

with open(infile, "r") as filelist:
    for line in filelist:
        filename, md5sum = line.strip().split()
        MATCHING = check_md5sum(filename, md5sum) and MATCHING

if not MATCHING:
    logging.critical("File contains mismatching md5sums")
    sys.exit(1)
