import json
import os
import sys
import argparse
import logging
import re

from collections import defaultdict

sys.path.insert(
    1, os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "lims",
        "stamlims_api"
))

from stamlims_api import rest

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
log = logging.getLogger(__name__)

rest.DEFAULT_ITEM_LIMIT = 10000
API = rest.setup_api({rest.RAISE_ON_ERROR_VAR: True})


lib_to_pool = defaultdict(list)

def parse_json(filename):
    with open(filename) as f:
        return json.loads(f.read())


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=parse_json)
    parser.add_argument("--output", default="sample_config.tsv")
    return parser


def group_data(processing_info) -> dict:
    """
    group_data tries to estimate what library pools each library belongs to
    Returns dict of tuple keys, values are a list of library numbers
    """
    output = defaultdict(list)
    for lib in processing_info['libraries']:
        lib_number = lib['library']
        key = (
                lib['barcode1']['reverse_sequence'],
                lib['lane'],
            )
        output[key].append(lib_number)

    return output

def populate_lib_to_pool():
    global lib_to_pool
    url_regex = re.compile("(\d+)")
    for pool in API.get_list_result(
            url_addition="library_pool/",
            item_limit=10000):
        name = pool['object_name']
        for lib_url in pool['libraries']:
            match = url_regex.search(lib_url)
            if not match:
                raise Exception("lib url %s didn't match" % lib_url)
            lib_id = int(match.group(1))
            lib_to_pool[lib_id].append(name)



# Ugly hack lol
# TODO: Add pool to processing_info endpoint, then we can remove this.
def get_pools_for_libs(lib_numbers) -> set:
    """ Returns dict of {pool: [lib_numbers]} """
    global lib_to_pool
    pools = set()

    numbers = ",".join(str(n) for n in lib_numbers)
    url_addition = "library/"
    for lib in API.get_list_result(
        url_addition=url_addition,
        item_limit=200,
        query_arguments={"number__in": numbers}
    ):
        id = int(lib['id'])
        assert len(lib_to_pool[id]) == 1, "Library LN%d should have exactly one pool" % lib['number']
        pool = lib_to_pool[id][0]
        pools.add(pool)

    return pools
        

def to_tsv(label, data):
    lines = ["name\tlane\tbarcode_index"]
    for (pool, index, lane, _numbers)  in sorted(data, key=lambda x:  (x[2],x[0])):
        lines.append(
                "\t".join([label + "_" + pool, str(lane), index])
        )
    return "\n".join(lines) + "\n"


# def create_upload_script(label, data):
#     
#     base = 'python3 "$STAMPIPES/scripts/lims/upload_data.py --attach_file_contenttype SequencingData.flowcelllane "
#     lines = ["#!/bin/bash"]
#     for (prefix, lane_ids) in data:
#         for num in numbers:
#             r1 = 
#             lines.append(
#                     base + " --attach_file_objectid %d --attach_file %s --attach-file-purpose  r1-fastq --attach-file-type fastq" % (num


def main():
    poptions = parser_setup().parse_args()
    label = poptions.data["flowcell"]["label"]
    grouped = group_data(poptions.data)

    populate_lib_to_pool()

    output_data = []
    for group_key, numbers in grouped.items():
        pools = get_pools_for_libs(numbers)
        pool_name = "_and_".join(sorted(pools))
        output_data.append( (pool_name, group_key[0], group_key[1], numbers) )


    tsv = to_tsv(label, output_data)
    with open(poptions.output, 'w') as f:
        f.write(tsv)

    # upload_data = []
    # for library in output_data["libraries"]:
    #     upload_data.append(library['id'], library['number'])
    # upload_script = create_upload_script(label, output_data)

if __name__ == "__main__":
    main()
