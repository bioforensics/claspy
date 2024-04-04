# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------
from claspy.db import CellosaurusDB
from claspy import Profile
from claspy.tests import data_file
from io import StringIO
import pytest


def test_search_report_sorting():
    db = CellosaurusDB.load(data_file("snu-db.json"))
    query = next(Profile.load(data_file("snu-query.csv")))
    report = db.search(query).full_report
    print(report.to_string())
    assert report.CellLine.to_list() == ["sample", "SNU-1033-1", "SNU-1033-3", "SNU-1033-2"]
    for score in report.Score[1:]:
        assert score == pytest.approx(1.0)
    assert report.SharedAlleles[1:].to_list() == [29, 27, 20]


def test_db_round_trip(tmp_path):
    db = CellosaurusDB.load()
    db1 = CellosaurusDB([profile for profile in db if "SK-HEP-1" in profile.identifier])
    db1.to_json(tmp_path / "db.json")
    db2 = CellosaurusDB.load(tmp_path / "db.json")
    assert len(db1) == len(db2)
    json1, json2 = StringIO(), StringIO()
    db1.to_json(json1)
    db2.to_json(json2)
    assert json1.getvalue() == json2.getvalue()
