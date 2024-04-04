# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .str_profile import Profile
from .result import ProfileResult, SearchResult
from importlib.resources import files
import json
from pathlib import Path
import re
import sys
from tqdm import tqdm
from urllib.request import urlretrieve


class CellosaurusDB(list):
    def search(
        self,
        query,
        algorithm="Tanabe",
        mode="intersect",
        amel=False,
        taxid=9606,
        minscore=0.0,
        maxhits=20,
    ):
        result = SearchResult(query, minscore=minscore, maxhits=maxhits)
        for reference in self:
            if taxid is not None and not reference.taxid_match(taxid):
                continue
            score, num_shared_alleles = Profile.score(
                query, reference, algorithm=algorithm, mode=mode, amel=amel
            )
            proresult = ProfileResult(query._meta["sample"], score, num_shared_alleles, reference)
            result.add_profile_result(proresult)
        return result

    @classmethod
    def load(cls, path=None):
        if path is None:
            path = cls.default_path()
        with open(path, "r") as instream:
            return cls.from_json(instream)

    @staticmethod
    def default_path():
        return files("claspy") / "cellosaurus.json"

    @classmethod
    def from_json(cls, instream):
        payload = json.load(instream)
        if not isinstance(payload, dict) and not isinstance(payload, list):
            raise ValueError(f"unexpected data type '{type(payload)}'")
        if isinstance(payload, dict):
            payload = [payload]
        records = cls()
        for profile in payload:
            metadata = profile["meta"]
            alleles = profile["alleles"]
            records.append(Profile(alleles, metadata))
        return records

    @classmethod
    def convert_from_download(cls, url=None):
        if url is None:
            url = "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt"
        path = files("claspy") / "cellosaurus.txt"
        with ProgressBar(unit="B", unit_scale=True, miniters=1, desc=Path(url).name) as pb:
            urlretrieve(url, path, reporthook=pb.update_to)
        return cls.convert_from_path(path)

    @classmethod
    def convert_from_path(cls, path=None):
        profiles = cls()
        with open(path, "r") as instream:
            parser = cls.parse_cellosaurus_records(instream)
            for n, profile in enumerate(parser):
                profiles.append(profile)
        print(f"[CellosaurusDB] parsed {n+1} distinct cell line STR profiles", file=sys.stderr)
        return profiles

    @staticmethod
    def parse_cellosaurus_records(instream):
        parser = CellosaurusDB.parse_cellosaurus_into_blocks(instream)
        for n, block in enumerate(parser):
            entry = CellosaurusEntry(block)
            for alleles, meta in entry.profiles:
                yield Profile(alleles, meta)
        print(f"[CellosaurusDB] parsed {n+1} database records", file=sys.stderr)

    @staticmethod
    def parse_cellosaurus_into_blocks(instream):
        block = list()
        for line in instream:
            if line.startswith("ID"):
                break
        block.append(line.strip())
        for line in instream:
            line = line.strip()
            if line == "//":
                yield block
                block = list()
            block.append(line)

    def to_json(self, output):
        if isinstance(output, str) or isinstance(output, Path):
            with open(output, "w") as outstream:
                json.dump([profile.payload for profile in self], outstream, indent=4)
        else:
            json.dump([profile.payload for profile in self], output, indent=4)


class CellosaurusEntry:
    ATTRIBUTES = {
        "ID": "identifier",
        "AC": "accession",
        "SY": "synonyms",
    }

    def __init__(self, data):
        self._data = data
        self.meta = dict()
        self.alleles = dict()
        for line in data:
            self.parse_meta(line)
            self.parse_sources(line)
            self.parse_alleles(line)

    def parse_meta(self, line):
        if line.startswith(("ID", "AC", "SY")):
            key, value = re.split(r"\s+", line, 1)
            assert key not in self.meta, key
            self.meta[self.ATTRIBUTES[key]] = value
        elif line.startswith("OX"):
            match = re.match(r"OX   NCBI_TaxID=(\d+); ! ([^\n]+)", line)
            if not match:
                raise ValueError(f"cannot parse species of origin: {line}")
            taxid, organism = match.groups()
            if "taxid" not in self.meta:
                self.meta["taxid"] = list()
                self.meta["organism"] = list()
            self.meta["taxid"].append(int(taxid))
            self.meta["organism"].append(organism)

    def parse_sources(self, line):
        if line.startswith("ST") and "Source" in line:
            match = re.match(r"ST   Source\(s\): ([^\n]+)", line)
            if not match:
                raise ValueError(f"could not parse sources: {line}")
            source_string = match.group(1)
            for source in source_string.split("; "):
                self.alleles[source] = dict()

    def parse_alleles(self, line):
        if line.startswith("ST") and "Source" not in line and "Not_detected" not in line:
            match = re.match(r"^ST   ([^:]+): ([\dXY,\. ]+)(.+)?", line)
            if not match:
                raise ValueError(f"could not parse STR profile data: {line}")
            marker, allele_str, sources = match.groups()
            if sources is None:
                for marker_alleles in self.alleles.values():
                    marker_alleles[marker] = allele_str.strip()
            else:
                sources = sources.replace("(", "").replace(")", "")
                for source in sources.split("; "):
                    if source not in self.alleles:
                        print(
                            "[CellosaurusDB] WARNING:",
                            f"Source '{source}' not defined for cell line {self.meta['identifier']}",
                            file=sys.stderr,
                        )
                    else:
                        self.alleles[source][marker] = allele_str.strip()

    @property
    def profiles(self):
        for source, marker_alleles in self.alleles.items():
            metadata = dict(self.meta)
            if len(metadata["taxid"]) == 1:
                metadata["taxid"] = metadata["taxid"][0]
                metadata["organism"] = metadata["organism"][0]
            metadata["source"] = source
            yield marker_alleles, metadata


class ProgressBar(tqdm):
    """Stolen shamelessly from https://stackoverflow.com/a/53877507/459780."""

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)
