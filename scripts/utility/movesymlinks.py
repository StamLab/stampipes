#!/usr/bin/env python

"""
Use to update symlinks with a target base path of X to the new target
base path of Y.
"""

import argparse
import logging
import os
import sys

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"


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
        "--fromdir",
        default=os.getcwd(),
        help="The directory we're changing all the symlinks in, defaults to current working directory",
    )
    parser.add_argument(
        "--olddir",
        required=True,
        help="The old directory we don't want symlinks going to",
    )
    parser.add_argument(
        "--newdir",
        required=True,
        help="The new base directory the symlinks should go to",
    )
    parser.add_argument(
        "--report", help="The outfile to write down all symlinks", default=None
    )
    parser.add_argument(
        "--move", dest="move", action="store_true", help="Actually perform the move"
    )
    return parser


class SymlinkMover(object):
    def __init__(self, fromdir, olddir, newdir, move=False, report=None):
        self.fromdir = fromdir
        self.olddir = olddir
        self.newdir = newdir
        self.domove = move
        self.alllinks = []
        self.movedlinks = []
        self.brokenlinks = []

    def detect(self, path):
        path = path.rstrip("/")

        if not os.path.islink(path):
            logging.debug("%s not a symlink" % path)
            return

        target_path = os.readlink(path)
        broken = False
        logging.debug("checking %s" % path)
        # Resolve relative symlinks
        if not os.path.isabs(target_path):
            target_path_absolute = os.path.join(os.path.dirname(path), target_path)
        else:
            target_path_absolute = target_path
        self.alllinks.append(path)
        if self.olddir in target_path:
            logging.debug("path %s target %s" % (path, target_path_absolute))
            self.movedlinks.append(path)
        if not os.path.exists(target_path_absolute):
            broken = True
            self.brokenlinks.append(path)

        if self.report:
            self.report.write(
                "%s\t%s\t%s\t%s\n"
                % (path, target_path, target_path_absolute, str(broken))
            )

    def move_link(self, linkpath):
        old_target_path = os.readlink(linkpath)
        new_target_path = old_target_path.replace(self.olddir, self.newdir)

        logging.info(
            "Moving %s pointer from %s to %s"
            % (linkpath, old_target_path, new_target_path)
        )
        try:
            if self.domove:
                os.unlink(linkpath)
                os.symlink(new_target_path, linkpath)
        except PermissionError:
            logging.error("Couldn't move %s, permission denied" % linkpath)

    def walk(self, directory):
        for root, dirs, files in os.walk(directory):
            if root.startswith("./.git"):
                # Ignore the .git directory.
                continue
            logging.debug("walking through directories for %s" % root)
            [self.detect(os.path.join(root, dirname)) for dirname in dirs]
            logging.debug("walking through files for %s" % root)
            [self.detect(os.path.join(root, filename)) for filename in files]

    def run(self, report=None):
        if report:
            self.report = open(report, "w")
        else:
            self.report = None

        logging.info("Detecting symlinks")
        self.walk(self.fromdir)
        logging.info("%d symlinks found in total" % len(self.alllinks))
        if self.brokenlinks:
            logging.info("broken symlink(s) found:")
            for link in self.brokenlinks:
                logging.info("\t%s" % link)
        if self.movedlinks:
            logging.info("%s symlinks to move" % len(self.movedlinks))
            logging.info("Symlink moves...")
            [self.move_link(link) for link in self.movedlinks]

        if self.report:
            self.report.close()


def main(args=sys.argv):
    parser = parser_setup()
    args = parser.parse_args()

    if args.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif args.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)

    mover = SymlinkMover(args.fromdir, args.olddir, args.newdir, args.move)
    mover.run(report=args.report)


if __name__ == "__main__":
    main()
