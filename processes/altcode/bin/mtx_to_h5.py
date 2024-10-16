#!/usr/bin/env python3

import argparse
import json
import logging

import scanpy as sc


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("mtx_directory", help="the directory containing the mtx files")
    parser.add_argument("output", help="the name of the output file")
    parser.add_argument("--metadata", help="A JSON-formatted file of metadata to add")
    parser.add_argument(
        "--compress", action="store_true", help="Compress output with gzip"
    )
    return parser


def lists_to_dicts(data):
    # Recursively converts lists to dicts
    # Required because of this issue https://github.com/scverse/anndata/issues/708
    if isinstance(data, list):
        return {f"_{idx}": lists_to_dicts(elem) for idx, elem in enumerate(data)}
    if isinstance(data, dict):
        for key in list(data.keys()):
            data[key] = lists_to_dicts(data[key])
    return data


def convert(input_dir, output_file, compress=False, metadata=None):
    data = sc.read_10x_mtx(input_dir, cache=False)
    if metadata is not None:
        # We store it two different ways because AnnData does not support lists-of-dicts
        data.uns["metadata_json"] = json.dumps(metadata)
        data.uns["metadata"] = lists_to_dicts(metadata)
    comp_method = "gzip" if compress else None
    data.write(filename=output_file, compression=comp_method)


def main():
    poptions = parser_setup().parse_args()
    if not poptions.output.endswith("h5ad"):
        logging.warning(
            "output file extension is not '.h5ad', some programs may fail to read it"
        )

    if poptions.metadata:
        with open(poptions.metadata) as m:
            metadata = json.load(m)
    else:
        metadata = None

    convert(
        poptions.mtx_directory, poptions.output, poptions.compress, metadata=metadata
    )


if __name__ == "__main__":
    main()
