# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
# This file is part of claspy: https://github.com/bioforensics/claspy
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------


valid_names = {
    9606: [
        "Amelogenin",
        "CSF1PO",
        "D10S1248",
        "D12S391",
        "D13S317",
        "D16S539",
        "D17S1301",
        "D18S51",
        "D19S433",
        "D1S1656",
        "D20S482",
        "D21S11",
        "D22S1045",
        "D2S1338",
        "D2S441",
        "D3S1358",
        "D4S2408",
        "D5S818",
        "D6S1043",
        "D7S820",
        "D8S1179",
        "D9S1122",
        "DXS10074",
        "DXS101",
        "DXS10103",
        "DXS10135",
        "DXS7132",
        "DXS7423",
        "DXS8378",
        "DYF387S1",
        "DYS19",
        "DYS385a-b",
        "DYS389I",
        "DYS389II",
        "DYS390",
        "DYS391",
        "DYS391",
        "DYS392",
        "DYS437",
        "DYS438",
        "DYS439",
        "DYS448",
        "DYS460",
        "DYS481",
        "DYS505",
        "DYS522",
        "DYS533",
        "DYS549",
        "DYS570",
        "DYS570",
        "DYS576",
        "DYS576",
        "DYS612",
        "DYS635",
        "DYS643",
        "F13A01",
        "F13B",
        "FESFPS",
        "FGA",
        "HPRTB",
        "LPL",
        "Penta C",
        "Penta D",
        "Penta E",
        "SE33",
        "TH01",
        "TPOX",
        "Y-GATA-H4",
        "vWA",
    ],
    10090: [
        "Mouse STR 1-1",
        "Mouse STR 1-2",
        "Mouse STR 2-1",
        "Mouse STR 3-2",
        "Mouse STR 4-2",
        "Mouse STR 5-5",
        "Mouse STR 6-4",
        "Mouse STR 6-7",
        "Mouse STR 7-1",
        "Mouse STR 8-1",
        "Mouse STR 9-2",
        "Mouse STR 11-2",
        "Mouse STR 12-1",
        "Mouse STR 13-1",
        "Mouse STR 15-3",
        "Mouse STR 17-2",
        "Mouse STR 18-3",
        "Mouse STR 19-2",
        "Mouse STR X-1",
    ],
    9615: [
        "Dog FHC2010",
        "Dog FHC2054",
        "Dog FHC2079",
        "Dog PEZ1",
        "Dog PEZ3",
        "Dog PEZ5",
        "Dog PEZ6",
        "Dog PEZ8",
        "Dog PEZ12",
        "Dog PEZ20",
    ],
}

species_by_taxid = {
    9606: "human",
    10090: "mouse",
    9615: "dog",
}


def validate_names(marker_names):
    """Validate marker names

    For each given marker name, determine the standardized form. Determine the species associated
    with this list of marker names. Raise an exception if any marker name cannot be validated, or
    if the list contains marker names from multiple species.
    """
    taxids = set()
    valid = dict()
    for name in marker_names:
        valid[name], taxid = standardize_name(name)
        taxids.add(taxid)
    if None in taxids:
        invalid = [name for name, valid_name in valid.items() if valid_name is None]
        invalid = ", ".join(invalid)
        raise ValueError(f"invalid marker name(s): {invalid}")
    if len(taxids) > 1:
        species = sorted([species_by_taxid[taxid] for taxid in taxids])
        species = ", ".join(species)
        message = f"list of marker names includes markers from different species: {species}"
        raise ValueError(message)
    taxid = taxids.pop()
    return valid, taxid


def standardize_name(name):
    candidate = name.replace(" ", "").lower()
    for taxid, species_names in valid_names.items():
        for species_name in species_names:
            species_candidate = species_name.replace(" ", "").lower()
            if candidate == species_candidate:
                return species_name, taxid
    return None, None
