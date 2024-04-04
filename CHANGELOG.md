# Change Log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Unreleased

### Changed
- Improvements to loading profile and database objects (!9)
- Database search is now restricted based on species inferred from markers in the query profile, not by user-specified species (!8)
- Summary report is displayed in terminal, full report to a CSV file (!11, !12)

### Fixed
- Added names of additional valid markers present in ForenSeq but not in Cellosaurus; includes four autosomal, seven X chromosome, and 21 Y chromosome STR markers (!8)
- Rank order of results with the same score but different numbers of shared alleles (!10)


## [0.0.2] 2023-05-25

### Fixed
- Divide by zero bug when query and reference have no shared alleles (!6)
- Marker name validation and standardization for human, mouse, and dog (!6)
- Rank order of results with the same score but different numbers of shared alleles (!6)
- Handling of string alleles, e.g. X and Y for Amelogenin (!7)
- Smart natural (not lexicographical) sorting of alleles for display (!7)


## [0.0.1] 2023-05-22

Initial release! Includes:

- `claspy_db` for downloading and formatting the Cellosaurus database
- `claspy` for searching a profile against Cellosaurus and reporting the best results
