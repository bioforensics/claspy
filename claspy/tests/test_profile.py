# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from claspy import Profile
from claspy.db import CellosaurusDB
from claspy.tests import data_file
import pytest


def test_profile_basic(capsys):
    alleles = {
        "CSF1PO": "13,14",
        "D5S818": "13",
        "D7S820": "8",
        "D13S317": "12",
        "FGA": "24",
        "TH01": "8",
        "TPOX": "11",
        "vWA": "16",
    }
    meta = {"sample": "sample1"}
    profile = Profile(alleles, meta)
    assert len(profile) == 9
    assert next(iter(profile)) == ("CSF1PO", "13")
    assert profile.taxid == 9606
    profile._meta["taxid"] = [9606, 10116]
    assert profile.taxid_match(9606) is True
    assert profile.taxid_match(10116) is True
    assert profile.taxid_match(10090) is False
    score, num_shared_alleles = Profile.score(profile, profile)
    assert score == pytest.approx(1.0)
    assert num_shared_alleles == 9
    print(profile)
    terminal = capsys.readouterr()
    observed = terminal.out
    expected = """
Sample,Marker,Allele1,Allele2
sample1,CSF1PO,13,14
sample1,D5S818,13,
sample1,D7S820,8,
sample1,D13S317,12,
sample1,FGA,24,
sample1,TH01,8,
sample1,TPOX,11,
sample1,vWA,16,
"""
    assert observed.strip() == expected.strip()


@pytest.mark.parametrize(
    "path,message",
    [
        ("query-wide.csv", r"found 11 allele columns, well above expected limit"),
        ("query-bad-allele-1.csv", r"expected column 'Allele1' missing"),
        ("query-bad-allele-2.csv", r"invalid table header 'AlleleTwo'"),
    ],
)
def test_load_failure_modes(path, message):
    with pytest.raises(ValueError, match=message):
        next(Profile.load(data_file(path)))


@pytest.mark.parametrize(
    "input,expected",
    [
        ("7", {"7"}),
        ("7,8", {"7", "8"}),
        ("7, 8", {"7", "8"}),
        ("10.2,13", {"10.2", "13"}),
        ("11.1", {"11.1"}),
    ],
)
def test_parse_allele_string(input, expected):
    assert Profile.parse_allele_string(input) == expected


@pytest.mark.parametrize(
    "input,expected",
    [
        ({"13", "7"}, "7,13"),
        ({"13"}, "13"),
        ({"19", "21", "9"}, "9,19,21"),
    ],
)
def test_allele_representation_is_sorted(input, expected):
    assert Profile.allele_repr(input) == expected


def test_allele_transform_failure_mode():
    with pytest.raises(ValueError, match=r"unexpected allele 'Z'"):
        Profile.allele_transform("Z")


@pytest.mark.parametrize(
    "algorithm,mode,amel,exp_alleles,exp_score",
    [
        ("Tanabe", "intersect", True, 21, 0.933333),
        ("Tanabe", "intersect", False, 19, 0.9268293),
        ("Tanabe", "query", False, 19, 0.904762),
        ("Tanabe", "reference", False, 19, 0.883721),
        ("query", "intersect", False, 19, 0.95),
        ("reference", "intersect", False, 19, 0.904762),
    ],
)
def test_score_basic(algorithm, mode, amel, exp_alleles, exp_score):
    profiles = CellosaurusDB.load(path=data_file("upci-scc-077-db.json"))
    assert len(profiles) == 2
    query, reference = profiles
    score, shared_alleles = Profile.score(
        query, reference, algorithm=algorithm, mode=mode, amel=amel
    )
    assert score == pytest.approx(exp_score)
    assert shared_alleles == exp_alleles


@pytest.mark.parametrize(
    "algorithm,mode,message",
    [
        ("Tanabe", "lizard", r"unsupported scoring mode 'lizard'"),
        ("AI", "intersect", r"unsupported scoring algorithm 'AI'"),
    ],
)
def test_score_failure_modes(algorithm, mode, message):
    query = next(Profile.load(data_file("mock-cvcl-1085.csv")))
    reference = next(Profile.load(data_file("db-cvcl-1085.csv")))
    with pytest.raises(ValueError, match=message):
        Profile.score(query, reference, algorithm=algorithm, mode=mode)


def test_parse_and_score():
    profiles = CellosaurusDB.load(path=data_file("examples.json"))
    assert len(profiles) == 2
    score, num_shared_alleles = Profile.score(profiles[0], profiles[1])
    assert score == pytest.approx(0.7)
    assert num_shared_alleles == 7


def test_score():
    query = next(Profile.load(data_file("mock-cvcl-1085.csv")))
    reference = next(Profile.load(data_file("db-cvcl-1085.csv")))
    score, num_shared_alleles = Profile.score(query, reference)
    assert score == pytest.approx(0.9677, abs=1e-4)
    assert num_shared_alleles == 30


@pytest.mark.parametrize(
    "allele_set,expected",
    [
        ({"9", "13"}, "9,13"),
        ({"3", "21"}, "3,21"),
        ({"7"}, "7"),
        ({"Y", "X"}, "X,Y"),
        ({"9.3", "13"}, "9.3,13"),
    ],
)
def test_allele_repr(allele_set, expected):
    assert Profile.allele_repr(allele_set) == expected
