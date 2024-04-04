# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from .db import CellosaurusDB
from .str_profile import Profile
from argparse import ArgumentParser
from claspy.db import CellosaurusDB
import claspy
import pandas as pd
import sys


def main(arglist=None):
    if arglist:
        arglist = map(str, arglist)
    args = get_parser().parse_args(arglist)
    db = CellosaurusDB.load(args.db)
    all_summaries = list()
    all_reports = list()
    for query in Profile.load(args.query):
        results = db.search(
            query,
            algorithm=args.algorithm,
            mode=args.mode,
            taxid=query.taxid,
            amel=args.amel,
            minscore=args.min_score,
            maxhits=args.max_hits,
        )
        all_summaries.append(results.summary)
        all_reports.append(results.full_report)
    pd.concat(all_summaries).to_markdown(sys.stdout, index=False, floatfmt=".3f")
    print("")
    if args.out:
        pd.concat(all_reports).to_csv(args.out, index=False)
        print(f"\nFull report written to {args.out}", file=sys.stderr)


def get_parser():
    parser = ArgumentParser(description="Claspy: cell line authentication with STRs in Python")
    parser.add_argument("query", help="query STR profile")
    parser.add_argument(
        "-v", "--version", action="version", version=f"Claspy v{claspy.__version__}"
    )
    parser.add_argument(
        "-d",
        "--db",
        metavar="PATH",
        default=CellosaurusDB.default_path(),
        help=f"path to Cellosaurus database; default is {CellosaurusDB.default_path()}",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        metavar="A",
        choices=("Tanabe", "reference", "query"),
        default="Tanabe",
        help="scoring algorithm; available options are Tanabe (2S/(Q+R)), query (S/Q), and reference (S/R); default is Tanabe",
    )
    parser.add_argument(
        "-m",
        "--mode",
        metavar="M",
        choices=("intersect", "reference", "query"),
        default="intersect",
        help="mode for handling missing data; available options are query (all query markers), reference (all reference markers), and intersect (only shared markers); default is intersect",
    )
    parser.add_argument(
        "-s",
        "--min-score",
        type=float,
        metavar="S",
        default=0.0,
        help="do not report candidate matches with a score < S; by default S=0 (filter disabled)",
    )
    parser.add_argument(
        "-x",
        "--max-hits",
        type=int,
        metavar="X",
        default=20,
        help="do not report more than X candidate matches; by default X=20; set X<=0 to disable this filter",
    )
    parser.add_argument(
        "--amel",
        action="store_true",
        help="include the Amelogenin marker, if present, in scoring calculations; by default it is excluded",
    )
    parser.add_argument(
        "-o",
        "--out",
        metavar="FILE",
        help="write a full report in CSV format to FILE; by default only a summary report is printed to the terminal",
    )
    return parser


def db_main(arglist=None):
    args = get_db_parser().parse_args(arglist)
    if args.path is None:
        records = CellosaurusDB.convert_from_download()
    else:
        records = CellosaurusDB.convert_from_path(args.path)
    records.to_json(args.dest)
    print(f"Database written to {args.dest}", file=sys.stderr)


def get_db_parser():
    parser = ArgumentParser(
        description="Retrieve, format, and install the Cellosaurus database for Claspy"
    )
    parser.add_argument(
        "-p",
        "--path",
        help="install the Cellosaurus database from local file PATH rather than a remote URL",
    )
    parser.add_argument(
        "-d",
        "--dest",
        metavar="PATH",
        default=CellosaurusDB.default_path(),
        help=f"destination for the Cellosaurus database in JSON format; by default PATH={CellosaurusDB.default_path()}",
    )
    return parser
