import random

from typing import List, Tuple
import pytest

from barcode_masks import get_barcode_masks


@pytest.mark.parametrize(
    "_name,read_len,index_len,lib_lengths,expected",
    [
        ("basic-paired", 75, 8, [(8, 8)], ["y75,i8,i8,y75"]),
        ("toolong-paired", 75, 10, [(8, 8)], ["y75,i8n2,i8n2,y75"]),
        ("tooshort-paired", 75, 8, [(10, 10)], ["y75,i8,i8,y75"]),
        ("mixed-paired", 75, 8, [(8, 8), (8, 0)], ["y75,i8,i8,y75", "y75,i8,n8,y75"]),
    ],
)
def test_expected_index_masks(_name, read_len, index_len, lib_lengths, expected):
    """Run some table-driven tests to make sure we get the right output"""
    data = make_processing_json(read_len, index_len, lib_lengths)
    actual = get_barcode_masks(data)
    assert set(actual) == set(expected)


def gen_barcode(length: int) -> str:
    """Generates a random string of letters of length 'length'"""
    return "".join([random.choice(["A", "C", "T", "G"]) for _ in range(length)])


def make_processing_json(
    read_len: int,
    index_len: int,
    lib_index_lengths: List[Tuple[int, int]],
) -> dict:
    """Creates a minimal "processing" data structure"""
    return {
        "flowcell": {
            "read_length": read_len,
            "index_length": index_len,
        },
        "libraries": [
            {
                "barcode1": {"sequence": gen_barcode(bc1)},
                "barcode2": {"sequence": gen_barcode(bc2)},
            }
            for (bc1, bc2) in lib_index_lengths
        ],
    }
