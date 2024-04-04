# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from claspy.markers import validate_names
import pytest


@pytest.mark.parametrize(
    "input_names,expected_names,expected_taxid",
    [
        (
            ("PENTAD", "AmeLoGenIN", "D21S11"),
            {"PENTAD": "Penta D", "AmeLoGenIN": "Amelogenin", "D21S11": "D21S11"},
            9606,
        ),
        (
            ("d8s1179", "PentaE", "Se33", "Tpox"),
            {"d8s1179": "D8S1179", "PentaE": "Penta E", "Se33": "SE33", "Tpox": "TPOX"},
            9606,
        ),
        (
            ("CSF1PO", "fgA", "D3S1358"),
            {"CSF1PO": "CSF1PO", "fgA": "FGA", "D3S1358": "D3S1358"},
            9606,
        ),
        (
            (" mousestr1-2 ", "MouseSTR8-1"),
            {" mousestr1-2 ": "Mouse STR 1-2", "MouseSTR8-1": "Mouse STR 8-1"},
            10090,
        ),
    ],
)
def test_validate_names_basic(input_names, expected_names, expected_taxid):
    observed_names, observed_taxid = validate_names(input_names)
    assert observed_names == expected_names
    assert observed_taxid == expected_taxid


def test_validate_names_invalid_marker():
    with pytest.raises(ValueError, match=r"invalid marker name\(s\): Penta G"):
        validate_names(("CSF1PO", "Penta G", "D2S1338"))


def test_validate_names_mixed_species():
    message = r"list of marker names includes markers from different species: dog, human"
    with pytest.raises(ValueError, match=message):
        validate_names(("vWA", "DogPEZ8"))
