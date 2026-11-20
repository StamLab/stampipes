"""Tests for make_samplesheets RunInfo parsing and OverrideCycles padding."""

from pathlib import Path

import pytest
from make_samplesheets import (
    BasesMask,
    IndexReadInfo,
    get_barcode_assignments,
    make_samplesheet_header,
    parse_runinfo,
    validate_runinfo_against_mask,
    write_samplesheets,
)


def write_runinfo(path: Path, reads: list[dict]) -> Path:
    read_xml = "\n".join(
        (
            f'      <Read Number="{i}" NumCycles="{r["cycles"]}" '
            f'IsIndexedRead="{r["indexed"]}" '
            f'IsReverseComplement="{r.get("rc", "N")}"/>'
        )
        for i, r in enumerate(reads, start=1)
    )
    path.write_text(
        f"""<?xml version="1.0"?>
<RunInfo Version="7">
  <Run Id="TEST" Number="1">
    <Reads>
{read_xml}
    </Reads>
  </Run>
</RunInfo>
"""
    )
    return path


def dual_index_runinfo(
    tmp_path: Path,
    index_cycles: int,
    index2_rc: bool = True,
    filename: str = "RunInfo.xml",
) -> Path:
    return write_runinfo(
        tmp_path / filename,
        [
            {"cycles": 211, "indexed": "N"},
            {"cycles": index_cycles, "indexed": "Y", "rc": "N"},
            {
                "cycles": index_cycles,
                "indexed": "Y",
                "rc": "Y" if index2_rc else "N",
            },
            {"cycles": 91, "indexed": "N"},
        ],
    )


def minimal_processing(
    bases_mask: str = "y211,i10,i10,y91",
    barcode_source_values: tuple[str, str] = ("AAAAAAAA", "CCCCCCCC"),
    reverse_values: tuple[str, str] = ("TTTTTTTT", "GGGGGGGG"),
    dual: bool = True,
) -> dict:
    bc1 = {
        "sequence": barcode_source_values[0],
        "reverse_sequence": reverse_values[0],
    }
    bc2 = None
    if dual:
        bc2 = {
            "sequence": barcode_source_values[1],
            "reverse_sequence": reverse_values[1],
        }
    return {
        "alignment_group": {"bases_mask": bases_mask},
        "libraries": [
            {
                "id": 1,
                "lane": 1,
                "samplesheet_name": "SampleA",
                "barcode1": bc1,
                "barcode2": bc2,
            }
        ],
    }


def test_parse_runinfo_index_flags(tmp_path):
    path = dual_index_runinfo(tmp_path, index_cycles=10, index2_rc=True)
    infos = parse_runinfo(path)
    assert infos == [
        IndexReadInfo(num_cycles=10, is_reverse_complement=False),
        IndexReadInfo(num_cycles=10, is_reverse_complement=True),
    ]


def test_parse_runinfo_missing_rc_defaults_false(tmp_path):
    path = write_runinfo(
        tmp_path / "RunInfo.xml",
        [
            {"cycles": 151, "indexed": "N"},
            {"cycles": 8, "indexed": "Y"},
            {"cycles": 8, "indexed": "Y"},
            {"cycles": 151, "indexed": "N"},
        ],
    )
    # Rewrite without IsReverseComplement attributes
    path.write_text(
        """<?xml version="1.0"?>
<RunInfo>
  <Run Id="OLD">
    <Reads>
      <Read Number="1" NumCycles="151" IsIndexedRead="N"/>
      <Read Number="2" NumCycles="8" IsIndexedRead="Y"/>
      <Read Number="3" NumCycles="8" IsIndexedRead="Y"/>
      <Read Number="4" NumCycles="151" IsIndexedRead="N"/>
    </Reads>
  </Run>
</RunInfo>
"""
    )
    infos = parse_runinfo(path)
    assert infos == [
        IndexReadInfo(num_cycles=8, is_reverse_complement=False),
        IndexReadInfo(num_cycles=8, is_reverse_complement=False),
    ]


def test_parse_runinfo_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_runinfo(tmp_path / "missing.xml")


def test_parse_runinfo_malformed(tmp_path):
    path = tmp_path / "RunInfo.xml"
    path.write_text("<RunInfo><broken>")
    with pytest.raises(ValueError, match="Malformed"):
        parse_runinfo(path)


def test_validate_runinfo_length_mismatch(tmp_path):
    path = dual_index_runinfo(tmp_path, index_cycles=10)
    infos = parse_runinfo(path)
    mask = BasesMask.parse("y211,i8,i8,y91")
    with pytest.raises(ValueError, match="length mismatch"):
        validate_runinfo_against_mask(mask, infos)


@pytest.mark.parametrize(
    "index_cycles,bc_lens,index2_rc,expected_mask,expected_override",
    [
        # Incident regression: 8bp barcodes on 10-cycle RC index2
        (10, (8, 8), True, "y211,i8n2,n2i8,y91", "Y211;I8N2;N2I8;Y91"),
        # Same lengths, no RC on index2 -> trailing pad on both
        (10, (8, 8), False, "y211,i8n2,i8n2,y91", "Y211;I8N2;I8N2;Y91"),
        # Full-length barcodes on 8-cycle RC index2 (reverse-invariant)
        (8, (8, 8), True, "y211,i8,i8,y91", "Y211;I8;I8;Y91"),
        # Full-length 10bp barcodes
        (10, (10, 10), True, "y211,i10,i10,y91", "Y211;I10;I10;Y91"),
        # Single-index style: barcode2 absent, pad entire index2
        (10, (8, 0), True, "y211,i8n2,n10,y91", "Y211;I8N2;N10;Y91"),
    ],
)
def test_override_cycles_padding(
    tmp_path, index_cycles, bc_lens, index2_rc, expected_mask, expected_override
):
    path = dual_index_runinfo(tmp_path, index_cycles, index2_rc=index2_rc)
    infos = parse_runinfo(path)
    root = BasesMask.parse(f"y211,i{index_cycles},i{index_cycles},y91")
    validate_runinfo_against_mask(root, infos)
    adjusted = root.adjust_for_barcode_lengths(*bc_lens, index_reads=infos)
    assert str(adjusted) == expected_mask
    assert adjusted.to_override_cycles() == expected_override


def test_adjust_requires_index_reads():
    mask = BasesMask.parse("y211,i10,i10,y91")
    with pytest.raises(ValueError, match="index_reads"):
        mask.adjust_for_barcode_lengths(8, 8)


@pytest.mark.parametrize(
    "barcode_source,expected_index,expected_index2",
    [
        ("sequence", "AAAAAAAA", "CCCCCCCC"),
        ("reverse_sequence", "TTTTTTTT", "GGGGGGGG"),
    ],
)
def test_barcode_source_selection(barcode_source, expected_index, expected_index2):
    data = minimal_processing()
    assigns = get_barcode_assignments(data, barcode_source)
    assert assigns[0]["barcode1"] == expected_index
    assert assigns[0]["barcode2"] == expected_index2


def test_inactive_index2_omits_mismatch_setting(tmp_path):
    path = dual_index_runinfo(tmp_path, index_cycles=10, index2_rc=True)
    infos = parse_runinfo(path)
    root = BasesMask.parse("y211,i10,i10,y91")
    adjusted = root.adjust_for_barcode_lengths(8, 0, index_reads=infos)
    assert adjusted.index_active_states == [True, False]
    header = make_samplesheet_header("Altius", "2026-01-01", adjusted, "1,1")
    assert "BarcodeMismatchesIndex1,1" in header
    assert "BarcodeMismatchesIndex2" not in header
    assert "OverrideCycles,Y211;I8N2;N10;Y91" in header


def test_write_samplesheets_incident_regression(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    runinfo = dual_index_runinfo(tmp_path, index_cycles=10, index2_rc=True)
    infos = parse_runinfo(runinfo)
    data = minimal_processing(
        bases_mask="y211,i10,i10,y91",
        barcode_source_values=("AACCGCGG", "CTAGCGCT"),
        reverse_values=("CCGCGGTT", "AGCGCTAG"),
    )
    assigns = get_barcode_assignments(data, "reverse_sequence")
    write_samplesheets(
        name="Altius",
        filename_template="SampleSheet.withmask.{mask}.csv",
        date="2026-09-14",
        root_mask=data["alignment_group"]["bases_mask"],
        assignments=assigns,
        index_reads=infos,
        mismatches="1,1",
    )
    out = tmp_path / "SampleSheet.withmask.y211,i8n2,n2i8,y91.csv"
    assert out.is_file()
    text = out.read_text()
    assert "OverrideCycles,Y211;I8N2;N2I8;Y91" in text
    assert "1,SampleA,CCGCGGTT,AGCGCTAG" in text
    # Old buggy filename must not be produced
    assert not (tmp_path / "SampleSheet.withmask.y211,i8n2,i8n2,y91.csv").exists()


@pytest.mark.parametrize(
    "index_cycles,umi,expected_override",
    [
        ((17, 8), 9, "Y211;I8U9;I8;Y91"),
        ((17, 10), 9, "Y211;I8U9;N2I8;Y91"),
        ((19, 10), 9, "Y211;I8U9N2;N2I8;Y91"),
        ((17, 10), 0, "Y211;I8N9;N2I8;Y91"),
    ],
)
def test_index1_umi_override_cycles(index_cycles, umi, expected_override):
    infos = [
        IndexReadInfo(index_cycles[0], False),
        IndexReadInfo(index_cycles[1], True),
    ]
    root = BasesMask.parse(f"y211,i{index_cycles[0]},i{index_cycles[1]},y91")
    adjusted = root.adjust_for_barcode_lengths(8, 8, index_reads=infos, umi1_len=umi)
    assert adjusted.to_override_cycles() == expected_override
    assert adjusted.index_active_states == [True, True]


def test_index1_umi_too_short():
    infos = [IndexReadInfo(10, False), IndexReadInfo(10, True)]
    with pytest.raises(ValueError, match="too short"):
        BasesMask.parse("y211,i10,i10,y91").adjust_for_barcode_lengths(
            8, 8, index_reads=infos, umi1_len=9
        )


@pytest.mark.parametrize(
    "barcode1_extra,expected_umi",
    [
        ({"umi_method": "xgen-udi-umi"}, 9),
        ({"category": "IDT xGen UDI-UMI Adapters", "umi_method": None}, 9),
        ({"category": "IDT xGen UDI-UMI Adapters", "umi_method": "xgen-udi-umi"}, 9),
        ({"category": "TruSeq", "umi_method": None}, 0),
        ({"category": "ThruPlex UMT", "umi_method": "thruplex"}, 0),
        ({}, 0),
    ],
)
def test_index1_umi_detection(barcode1_extra, expected_umi):
    data = minimal_processing()
    data["libraries"][0]["barcode1"].update(barcode1_extra)
    assert get_barcode_assignments(data, "sequence")[0]["umi1"] == expected_umi


def test_write_samplesheets_splits_umi_libraries(tmp_path, monkeypatch):
    """UMI and non-UMI libraries with equal barcode lengths get separate sheets."""
    monkeypatch.chdir(tmp_path)
    runinfo = write_runinfo(
        tmp_path / "RunInfo.xml",
        [
            {"cycles": 211, "indexed": "N"},
            {"cycles": 17, "indexed": "Y", "rc": "N"},
            {"cycles": 10, "indexed": "Y", "rc": "Y"},
            {"cycles": 91, "indexed": "N"},
        ],
    )
    data = minimal_processing(bases_mask="y211,i17,i10,y91")
    umi_lib = dict(data["libraries"][0], id=2, samplesheet_name="SampleUMI")
    umi_lib["barcode1"] = {
        "sequence": "ACGATCAG",
        "reverse_sequence": "CTGATCGT",
        "umi_method": "xgen-udi-umi",
    }
    umi_lib["barcode2"] = {"sequence": "GCGCATAT", "reverse_sequence": "ATATGCGC"}
    data["libraries"].append(umi_lib)

    write_samplesheets(
        name="Altius",
        filename_template="SampleSheet.withmask.{mask}.csv",
        date="2026-11-20",
        root_mask=data["alignment_group"]["bases_mask"],
        assignments=get_barcode_assignments(data, "reverse_sequence"),
        index_reads=parse_runinfo(runinfo),
        mismatches="1,1",
    )
    umi_sheet = (tmp_path / "SampleSheet.withmask.y211,i8u9,n2i8,y91.csv").read_text()
    plain_sheet = (tmp_path / "SampleSheet.withmask.y211,i8n9,n2i8,y91.csv").read_text()
    assert "OverrideCycles,Y211;I8U9;N2I8;Y91" in umi_sheet
    assert "1,SampleUMI,CTGATCGT,ATATGCGC" in umi_sheet
    assert "SampleA" not in umi_sheet
    assert "1,SampleA,TTTTTTTT,GGGGGGGG" in plain_sheet
    assert "SampleUMI" not in plain_sheet
