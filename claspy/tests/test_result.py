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
from claspy.str_profile import Profile
from claspy.tests import data_file
from io import StringIO
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def skhep_result():
    db = CellosaurusDB.load(data_file("skhep1-db.json"))
    query = next(Profile.load(data_file("mock-sk-hep-1.csv")))
    return db.search(query, maxhits=5, minscore=0.9)


@pytest.fixture(scope="session")
def skhep_result_multi_sample():
    db = CellosaurusDB.load(data_file("skhep1-db.json"))
    all_results = list()
    for query in Profile.load(data_file("mock-sk-hep-1-2samples.csv")):
        all_results.append(db.search(query, maxhits=3, minscore=0.9))
    return all_results


def test_search_result_basic(skhep_result):
    assert skhep_result.maxhits == 5
    assert skhep_result.minscore == pytest.approx(0.9)
    assert len(skhep_result.results_by_cell_line) == 11


def test_search_result_summary(skhep_result):
    observed = skhep_result.summary
    assert len(observed) == 5
    exp_data = StringIO(
        """
Sample,CellLine,Score,SharedAlleles,Source
mock,SK-HEP-1,0.981818,27,PubMed=25877200
mock,SK-HEP-1-Cas9-727,0.980392,25,CCRID
mock,SK-HEP-1-Cas9-726,0.961538,25,CCRID
mock,SK-HEP-1-Cas9-728,0.961538,25,CCRID
mock,SK-HEP-1-Cas9-729,0.961538,25,CCRID"""
    )
    expected = pd.read_csv(exp_data)
    pd.testing.assert_frame_equal(observed, expected, check_exact=False, rtol=1e-6)


def test_search_result_full_report(skhep_result):
    observed = skhep_result.full_report
    assert len(observed) == 7
    observed = observed.to_csv(sep=";", index=False)
    expected = """
Sample;CellLine;Status;Score;SharedAlleles;Source;Amelogenin;CSF1PO;D12S391;D13S317;D16S539;D18S51;D19S433;D21S11;D2S1338;D3S1358;D5S818;D6S1043;D7S820;D8S1179;FGA;Penta D;Penta E;TH01;TPOX;vWA
mock;mock;query;;;;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;14;17;13,14;;7,9;9;14,17
mock;SK-HEP-1;best;0.9818181818181818;27;PubMed=25877200;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;13;7,9;9;14,17
mock;SK-HEP-1;worst;0.9818181818181818;27;ATCC;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;;7,9;9;14,17
mock;SK-HEP-1-Cas9-727;only;0.9803921568627451;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock;SK-HEP-1-Cas9-726;only;0.9615384615384616;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31,32;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock;SK-HEP-1-Cas9-728;only;0.9615384615384616;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31,32;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock;SK-HEP-1-Cas9-729;only;0.9615384615384616;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31,32;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17"""
    assert observed.strip() == expected.strip()


def test_search_result_summary_multisamples(skhep_result_multi_sample):
    observed = pd.concat(
        [result.summary for result in skhep_result_multi_sample], ignore_index=True
    )
    assert len(observed) == 6
    exp_data = StringIO(
        """
Sample,CellLine,Score,SharedAlleles,Source
mock_1,SK-HEP-1,0.981818,27,PubMed=25877200
mock_1,SK-HEP-1-Cas9-727,0.980392,25,CCRID
mock_1,SK-HEP-1-Cas9-726,0.961538,25,CCRID
mock_2,SK-HEP-1,0.981818,27,PubMed=25877200
mock_2,SK-HEP-1-Cas9-727,0.980392,25,CCRID
mock_2,SK-HEP-1-Cas9-726,0.961538,25,CCRID"""
    )
    expected = pd.read_csv(exp_data)
    pd.testing.assert_frame_equal(observed, expected, check_exact=False, rtol=1e-6)


def test_search_result_full_report_multisamples(skhep_result_multi_sample):
    observed = pd.concat([result.full_report for result in skhep_result_multi_sample])
    assert len(observed) == 10
    observed = observed.to_csv(sep=";", index=False)
    expected = """
Sample;CellLine;Status;Score;SharedAlleles;Source;Amelogenin;CSF1PO;D12S391;D13S317;D16S539;D18S51;D19S433;D21S11;D2S1338;D3S1358;D5S818;D6S1043;D7S820;D8S1179;FGA;Penta D;Penta E;TH01;TPOX;vWA
mock_1;mock_1;query;;;;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;14;17;13,14;;7,9;9;14,17
mock_1;SK-HEP-1;best;0.9818181818181818;27;PubMed=25877200;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;13;7,9;9;14,17
mock_1;SK-HEP-1;worst;0.9818181818181818;27;ATCC;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;;7,9;9;14,17
mock_1;SK-HEP-1-Cas9-727;only;0.9803921568627451;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock_1;SK-HEP-1-Cas9-726;only;0.9615384615384616;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31,32;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock_2;mock_2;query;;;;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;14;17;13,14;;7,9;9;14,17
mock_2;SK-HEP-1;best;0.9818181818181818;27;PubMed=25877200;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;13;7,9;9;14,17
mock_2;SK-HEP-1;worst;0.9818181818181818;27;ATCC;X;11,12;18;8,12;12;13,15;;29,31;20,23;16;10,13;11;8,11;13,14;17;13,14;;7,9;9;14,17
mock_2;SK-HEP-1-Cas9-727;only;0.9803921568627451;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
mock_2;SK-HEP-1-Cas9-726;only;0.9615384615384616;25;CCRID;X;11,12;18;8,12;12;13,15;12,15.2;29,31,32;20,23;16;10,13;11;8,11;13,14;17;;13,21;7,9;9;14,17
"""

    assert observed.strip() == expected.strip()
