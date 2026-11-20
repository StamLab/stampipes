"""Tests for index-1 UMI cycles and index-2 padding order in make_samplesheets."""

import pytest
from make_samplesheets import BasesMask, get_barcode_assignments, write_samplesheets


@pytest.mark.parametrize(
    "mask,umi,pad_first,expected",
    [
        # xGen UDI-UMI on NovaSeq X, 17/8 cycles
        ("y211,i17,i8,y91", 9, True, "Y211;I8U9;I8;Y91"),
        # 17/10 cycles (shared with 10bp libraries): UMI kept, i5 padded first
        ("y211,i17,i10,y91", 9, True, "Y211;I8U9;N2I8;Y91"),
        # Spare cycles after the UMI are skipped
        ("y211,i19,i10,y91", 9, True, "Y211;I8U9N2;N2I8;Y91"),
        # Non-UMI library on the same run
        ("y211,i17,i10,y91", 0, True, "Y211;I8N9;N2I8;Y91"),
        # Old behaviour unchanged without the new options
        ("y211,i17,i10,y91", 0, False, "Y211;I8N9;I8N2;Y91"),
        ("y151,i8,i8,y151", 0, True, "Y151;I8;I8;Y151"),
    ],
)
def test_override_cycles(mask, umi, pad_first, expected):
    adjusted = BasesMask.parse(mask).adjust_for_barcode_lengths(
        8, 8, umi1_len=umi, index2_pad_first=pad_first
    )
    assert adjusted.to_override_cycles() == expected
    assert adjusted.index_active_states == [True, True]


def test_index2_pad_first_absent_barcode2():
    adjusted = BasesMask.parse("y211,i10,i10,y91").adjust_for_barcode_lengths(
        8, 0, index2_pad_first=True
    )
    assert adjusted.to_override_cycles() == "Y211;I8N2;N10;Y91"
    assert adjusted.index_active_states == [True, False]


def test_index1_umi_too_short():
    with pytest.raises(ValueError, match="too short"):
        BasesMask.parse("y211,i10,i10,y91").adjust_for_barcode_lengths(8, 8, umi1_len=9)


def _lib(lib_id, name, bc1_extra=None):
    bc1 = {"sequence": "AAAAAAAA", "reverse_sequence": "TTTTTTTT"}
    bc1.update(bc1_extra or {})
    return {
        "id": lib_id,
        "lane": 1,
        "samplesheet_name": name,
        "barcode1": bc1,
        "barcode2": {"sequence": "CCCCCCCC", "reverse_sequence": "GGGGGGGG"},
    }


@pytest.mark.parametrize(
    "bc1_extra,expected_umi",
    [
        ({"umi_method": "xgen-udi-umi"}, 9),
        ({"category": "IDT xGen UDI-UMI Adapters", "umi_method": None}, 9),
        ({"category": "IDT xGen UDI-UMI Adapters", "umi_method": "xgen-udi-umi"}, 9),
        ({"category": "TruSeq", "umi_method": None}, 0),
        ({"category": "ThruPlex UMT", "umi_method": "thruplex"}, 0),
        ({}, 0),
    ],
)
def test_umi_detection(bc1_extra, expected_umi):
    data = {"libraries": [_lib(1, "S", bc1_extra)]}
    assert get_barcode_assignments(data, True, True)[0]["umi1"] == expected_umi


def test_write_samplesheets_splits_umi_libraries(tmp_path, monkeypatch):
    """UMI and non-UMI libraries with equal barcode lengths get separate sheets."""
    monkeypatch.chdir(tmp_path)
    data = {
        "libraries": [
            _lib(1, "SamplePlain", {"category": "TruSeq"}),
            _lib(2, "SampleUMI", {"category": "IDT xGen UDI-UMI Adapters"}),
        ]
    }
    write_samplesheets(
        name="Altius",
        filename_template="SampleSheet.withmask.{mask}.csv",
        date="2026-11-20",
        root_mask="y211,i17,i10,y91",
        assignments=get_barcode_assignments(data, True, True),
        mismatches="1,1",
        index2_pad_first=True,
    )
    umi = (tmp_path / "SampleSheet.withmask.y211,i8u9,n2i8,y91.csv").read_text()
    plain = (tmp_path / "SampleSheet.withmask.y211,i8n9,n2i8,y91.csv").read_text()
    assert "OverrideCycles,Y211;I8U9;N2I8;Y91" in umi
    assert "1,SampleUMI,TTTTTTTT,GGGGGGGG" in umi
    assert "SamplePlain" not in umi
    assert "OverrideCycles,Y211;I8N9;N2I8;Y91" in plain
    assert "SampleUMI" not in plain
