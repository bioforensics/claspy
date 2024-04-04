# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .markers import validate_names
import pandas as pd
import re


class Profile:
    """Class for handling STR profiles

    Includes methods for loading, saving, and scoring genetic profiles based on short tandem repeat
    (STR) markers. Each profile contains a set of STR alleles. In mammalian cell lines, these are
    expected to be diploid; that is, there should be at most observed two alleles for each marker.
    """

    def __init__(self, alleles, meta):
        self._meta = meta
        valid_names, taxid = validate_names(alleles.keys())
        self.taxid = taxid
        self._alleles = dict()
        for marker, marker_alleles in alleles.items():
            marker = valid_names[marker]
            marker_alleles = Profile.parse_allele_string(marker_alleles)
            self._alleles[marker] = marker_alleles

    @classmethod
    def load(cls, path):
        types = {f"Allele{i+1}": str for i in range(10)}
        data = pd.read_csv(path, sep=None, engine="python", dtype=types)
        for column in ("Sample", "Marker", "Allele1"):
            if column not in data.columns:
                raise ValueError(f"expected column '{column}' missing")
        for sample_name, sample_data in data.groupby("Sample"):
            numalleles = Profile.num_alleles_from_table(sample_data)
            metadata = {"sample": sample_data.Sample.iloc[0]}
            alleles = dict()
            for i, row in sample_data.iterrows():
                marker_alleles = list()
                for n in range(numalleles):
                    allele = row[f"Allele{n+1}"]
                    if not pd.isna(allele):
                        marker_alleles.append(allele)
                alleles[row.Marker] = ",".join(sorted(marker_alleles))
            yield Profile(alleles, metadata)

    @staticmethod
    def num_alleles_from_table(table):
        count = 1
        for column in table:
            if column.startswith("Allele"):
                try:
                    number = int(column[6:])
                except ValueError as verr:
                    raise ValueError(f"invalid table header '{column}'") from verr
                if number > count:
                    count = number
        if count > 10:
            raise ValueError(f"found {count} allele columns, well above expected limit")
        return count

    @staticmethod
    def parse_allele_string(alleles):
        return set(alleles.replace(" ", "").split(","))

    @staticmethod
    def allele_repr(allele_set):
        alleles = [Profile.allele_transform(a) for a in allele_set]
        alleles = sorted(alleles)
        alleles = [str(a) for a in alleles]
        return ",".join(alleles)

    @staticmethod
    def allele_transform(allele):
        if "." in allele:
            return float(allele)
        elif re.match(r"^\d+$", allele):
            return int(allele)
        elif allele in ("X", "Y"):
            return allele
        else:
            raise ValueError(f"unexpected allele '{allele}'")

    @property
    def table(self):
        sample = self._meta["sample"] if "sample" in self._meta else "sample"
        alleles = list()
        for marker, marker_alleles in self._alleles.items():
            allele_repr = Profile.allele_repr(marker_alleles)
            sorted_alleles = allele_repr.split(",")
            entry = [sample, marker, *sorted_alleles]
            while len(entry) < self.max_num_alleles + 2:
                entry.append(None)
            alleles.append(entry)
        colnames = ["Sample", "Marker"] + [f"Allele{i+1}" for i in range(self.max_num_alleles)]
        return pd.DataFrame(alleles, columns=colnames)

    @property
    def max_num_alleles(self):
        return max([len(allele_set) for allele_set in self._alleles.values()])

    def marker_alleles(self, markers):
        for marker in markers:
            if marker not in self._alleles:
                yield pd.NA
            else:
                yield Profile.allele_repr(self._alleles[marker])

    def __iter__(self):
        for marker, alleles in self._alleles.items():
            for allele in sorted(alleles):
                yield marker, allele

    def __len__(self):
        return len([allele for allele in self])

    @staticmethod
    def score(query, reference, algorithm="Tanabe", mode="intersect", amel=False):
        """Compute a similarity score between two profiles

        The score is based on the number of alleles shared between the query profile and the
        reference profile. Three scoring algorithms are implemented as described below: "Tanabe" is
        the default (Q=# query alleles, R=# reference alleles, S=# shared alleles).

        - "Tanabe": 2S / (Q+R)
        - "query": S / Q
        - "reference: R / Q

        There are also three modes for handling missing allele data in one or both profiles: the
        "intersect" mode is used by default.

        - "intersect": consider alleles only at markers present in both profiles
        - "query": consider alleles for markers present in the query profile, even if missing from
          the reference profile
        - "reference": consider alleles for markers present in the reference profile, even if
          missing from the query profile

        The Amelogenin marker (amel) is used for sex determination and is typically excluded from
        scoring. Set `amel=True` to include.

        The Tanabe algorithm is described in doi:10.11418/jtca1981.18.4_329, while the query and
        reference algorithms are described in doi:10.1073/pnas.121616198
        """
        if algorithm not in ("Tanabe", "query", "reference"):
            raise ValueError(f"unsupported scoring algorithm '{algorithm}'")
        if mode not in ("intersect", "query", "reference"):
            raise ValueError(f"unsupported scoring mode '{mode}'")
        markers = Profile.markers_for_scoring(query, reference, mode=mode, amel=amel)
        query_alleles = len(query.alleles(markers=markers))
        refr_alleles = len(reference.alleles(markers=markers))
        shared_alleles = len(query.alleles(markers=markers) & reference.alleles(markers=markers))
        score = 0.0
        if shared_alleles > 0:
            if algorithm == "Tanabe":
                score = (2 * shared_alleles) / (query_alleles + refr_alleles)
            elif algorithm == "query":
                score = shared_alleles / query_alleles
            else:
                score = shared_alleles / refr_alleles
        return score, shared_alleles

    @staticmethod
    def markers_for_scoring(query, reference, mode="intersect", amel=False):
        if mode == "intersect":
            markers = query.markers & reference.markers
        elif mode == "query":
            markers = query.markers
        else:
            markers = reference.markers
        if not amel:
            markers = [m for m in markers if m != "Amelogenin"]
        return markers

    @property
    def markers(self):
        return set(self._alleles.keys())

    def alleles(self, markers=None):
        if markers is None:
            markers = self.markers
        return {(marker, allele) for marker, allele in self if marker in markers}

    def __str__(self):
        return self.table.to_csv(index=False)

    def taxid_match(self, taxid):
        if isinstance(self._meta["taxid"], list):
            for testid in self._meta["taxid"]:
                if int(testid) == int(taxid):
                    return True
            return False
        else:
            return int(self._meta["taxid"]) == int(taxid)

    @property
    def identifier(self):
        return self._meta["identifier"]

    @property
    def source(self):
        return self._meta["source"]

    @property
    def payload(self):
        return {"meta": self._meta, "alleles": self.allele_dict}

    @property
    def allele_dict(self):
        return {marker: ",".join(sorted(alleles)) for marker, alleles in self._alleles.items()}

    @property
    def slug(self):
        return self.identifier, self._meta["accession"], self.source

    def __lt__(self, other):
        return self.slug < other.slug

    @property
    def payload(self):
        return {"meta": self._meta, "alleles": self.allele_dict}
