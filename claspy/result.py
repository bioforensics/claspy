# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from collections import defaultdict, namedtuple
import pandas as pd


class SearchResult:
    """Result for database search of a single query profile

    The SearchResult includes a score for every profile of the relevant species in the database.
    Distinct profiles for the same cell line are stored in a single CellLineResult object,
    accessible by that cell line's identifier.
    """

    def __init__(self, query, minscore=0.0, maxhits=20):
        self.query = query
        self.minscore = minscore
        self.maxhits = maxhits
        self.results_by_cell_line = defaultdict(CellLineResult)

    def add_profile_result(self, result):
        self.results_by_cell_line[result.reference.identifier].append(result)

    @property
    def summary(self):
        colnames = ["Sample", "CellLine", "Score", "SharedAlleles", "Source"]
        summary = pd.DataFrame([result.summary for result in self], columns=colnames)
        return summary

    @property
    def full_report(self):
        entries = list()
        markers = self.all_markers
        entry = (
            self.query._meta["sample"],
            self.query._meta["sample"],
            "query",
            pd.NA,
            pd.NA,
            pd.NA,
            *self.query.marker_alleles(markers),
        )
        entries.append(entry)
        for result in self:
            for entry in result.full_report(markers):
                entries.append(entry)
        colnames = ["Sample", "CellLine", "Status", "Score", "SharedAlleles", "Source"] + markers
        return pd.DataFrame(entries, columns=colnames)

    def __iter__(self):
        for n, identifier in enumerate(self.ids_by_score):
            if self.maxhits > 0 and n >= self.maxhits:
                return
            result = self.results_by_cell_line[identifier]
            if result.top_score < self.minscore:
                return
            yield result

    @property
    def ids_by_score(self):
        sorted_results = sorted(
            self.results_by_cell_line.values(),
            key=lambda result: (result.top_score, result.top_score_shared_alleles),
            reverse=True,
        )
        for result in sorted_results:
            yield result.identifier

    @property
    def all_markers(self):
        """Determine all markers to report

        This includes any marker for which allele data is present in the query or at least one of
        the database profiles to be included in the final full report.
        """
        markers = set()
        for marker, allele in self.query.alleles():
            markers.add(marker)
        for result in self:
            for subresult in result:
                for marker, allele in subresult.reference.alleles():
                    markers.add(marker)
        return sorted(markers)


class CellLineResult(list):
    """A list of query search scores and database profiles from the same cell line

    This class is a list of ProfileResult objects, and essentially provides some convenience
    functions for handling one or more scored profiles for a cell line from a database search.
    """

    @property
    def top_score(self):
        return max([single_result.score for single_result in self])

    @property
    def top_score_shared_alleles(self):
        return max(
            [
                single_result.shared_alleles
                for single_result in self
                if single_result.score == self.top_score
            ]
        )

    @property
    def identifier(self):
        ids = [single_result.reference.identifier for single_result in self]
        assert len(set(ids)) == 1
        return ids[0]

    @property
    def sample(self):
        samples = [single_result.sample for single_result in self]
        assert len(set(samples)) == 1
        return samples[0]

    @property
    def summary(self):
        results = sorted(self, reverse=True)
        best = results[0]
        best = (best.score, best.shared_alleles, best.reference.source)
        return self.sample, self.identifier, *best

    def full_report(self, markers):
        results = sorted(self, reverse=True)
        best = results[0]
        status = "best" if len(results) > 1 else "only"
        yield self.sample, self.identifier, status, *best.full_report(markers)
        if len(results) > 1:
            worst = results[-1]
            yield self.sample, self.identifier, "worst", *worst.full_report(markers)


class ProfileResult(namedtuple("ProfileResult", "sample score shared_alleles reference")):
    """Score from comparing a query profile to a single database reference profile"""

    @property
    def summary(self):
        return self.score, self.shared_alleles, self.reference.source

    def full_report(self, markers):
        return (*self.summary, *self.reference.marker_alleles(markers))
