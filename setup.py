# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from setuptools import setup
import versioneer


with open("README.md", "r") as infile:
    longdesc = infile.read()

setup(
    name="claspy",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Clapsy: cell line authentication with STRs in Python",
    long_description=longdesc,
    long_description_content_type="text/markdown",
    url="https://github.com/bioforensics/claspy",
    packages=["claspy", "claspy.tests"],
    package_data={"claspy": ["claspy/tests/data/*"]},
    include_package_data=True,
    install_requires=[
        "black==24.3",
        "pandas>=2.0",
        "pytest>=6.0",
        "pytest-cov>=3.0",
        "tabulate>=0.9",
        "tqdm>=3.0",
    ],
    entry_points={"console_scripts": ["claspy = claspy:main", "claspy_db = claspy:db_main"]},
    classifiers=[
        "Environment :: Console",
        "Framework :: IPython",
        "Framework :: Jupyter",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=True,
)
