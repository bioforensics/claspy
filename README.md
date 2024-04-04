# Claspy: cell line authentication with STRs in Python

Documentation for Claspy is pending.
In the mean time, see the following hints.

```
claspy_db  # Run one time to install Cellosaurus database
claspy query.csv  # Run to find closest profile to the query in the database
```

STR profiles should be in tabular/CSV format and look something like this.

```csv
Sample,Marker,Allele1,Allele2
sample1,CSF1PO,12,13
sample1,D13S317,12,
sample1,D16S539,9,11
sample1,D18S51,12,15
sample1,D19S433,13,15
sample1,D21S11,29,32.2
sample1,D2S1338,20,23
sample1,D3S1358,16,17
sample1,D5S818,10,11
sample1,D7S820,10,11
sample1,D8S1179,13,15
sample1,FGA,18,24
sample1,Penta D,9,
sample1,Penta E,17,
sample1,TH01,9,9.3
sample1,TPOX,8,
sample1,vWA,18,19
```
