#!/usr/bin/env python3

import argparse
import logging

import scanpy as sc


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument("mtx_directory", help="the directory containing the mtx files")
    parser.add_argument("output", help="the name of the output file")
    parser.add_argument("--compress", action="store_true",
                        help="Compress output with gzip")
    return parser


def convert(input_dir, output_file, compress=False):
    data = sc.read_10x_mtx(input_dir, cache=False)
    comp_method = "gzip" if compress else None
    data.write(filename=output_file, compression=comp_method)

def main():
    poptions = parser_setup().parse_args()
    if not poptions.output.endswith("h5ad"):
        logging.warning(
            "output file extension is not '.h5ad', some programs may fail to read it"
        )
    convert(poptions.mtx_directory, poptions.output, poptions.compress)

if __name__ == "__main__":
    main()
