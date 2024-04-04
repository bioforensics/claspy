# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import claspy
from claspy.tests import data_file
import pytest


def test_search_report_sorting(tmp_path):
    report = tmp_path / "report.csv"
    arglist = [data_file("mock-cvcl-1085.csv"), "--out", report]
    claspy.cli.main(arglist=arglist)
    assert report.is_file()
    with open(report, "r") as fh1, open(data_file("report-cvcl-1085.csv"), "r") as fh2:
        observed = fh1.read().strip()
        expected = fh2.read().strip()
        assert observed == expected
