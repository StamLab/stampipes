import os, sys, logging, re
import requests
import json
import fileinput
import argparse

token = None
headers = None
lane_tags = None
flowcell_lane_cache = dict()
flowcell_contenttype = None
base_api_url = None

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

script_options = {
    "basedir": os.getcwd(),
    "quiet": False,
    "debug": False,
    "content_type": None,
    "object_id": None,
    "old_tag": None,
    "new_tag": None,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.")
    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required.")

    parser.add_argument("--content_type", dest="content_type",
        help="Name of the contenttype.")
    parser.add_argument("--object", dest="object_id", type=int,
        help="Object ID.")

    parser.add_argument("-a", "--add_tag", dest="new_tag",
        help="The new tag slug.")
    parser.add_argument("-r", "--remove_tag", dest="old_tag",
        help="The old tag slug.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

class TagChange(object):

    def __init__(self, api_url, token):
       self.api_url = api_url
       self.token = token
       self.headers = {'Authorization': "Token %s" % token}
       self.contenttypes = {}

    def get_contenttype(self, contenttype):
        if not contenttype in self.contenttypes:
            contenttype_url = '%s/content_type/?model=%s' % (self.api_url, contenttype)
            contenttype_results = requests.get(contenttype_url, headers = self.headers).json()
            self.contenttypes[contenttype] = contenttype_results['results'][0]
        return self.contenttypes[contenttype]

    def get_tag(self, slug):

        if not slug:
            return None

        exists = requests.get("%s/tag/?slug=%s" % (self.api_url, slug), headers = self.headers)
        tag = None
        if exists.ok:
            results = exists.json()
            if results['count'] > 0:
                return results['results'][0]
            else:
                print "Tag %s not found" % slug
                return None
        else:
            print "Error finding tag %s through API" % slug
            return None

    def change_tag(self, contenttype, object_id, old_tag, new_tag):

        contenttype_id = self.get_contenttype(contenttype)["id"]

        current = requests.get("%s/tagged_object/?content_type=%d&object_id=%d&tag__slug=%s" % (
            self.api_url, contenttype_id, object_id, old_tag), headers = self.headers).json()

        if current['count'] == 0:
            return False

        taggedobject = current['results'][0]
        new_tag = self.get_tag(new_tag)
        taggedobject["tag"] = new_tag["url"]

        result = requests.put(taggedobject['url'], headers = self.headers, data = taggedobject)

        if not result.ok:
            return False

        return True

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
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
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

    tagchange = TagChange(api_url, token)
    result = tagchange.change_tag(poptions.content_type, poptions.object_id, poptions.old_tag, poptions.new_tag)

    if not result:
        sys.stderr.write("Tag not changed")
        sys.exit(1)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
