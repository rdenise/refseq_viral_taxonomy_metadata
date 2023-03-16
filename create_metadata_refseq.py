#!/usr/bin/env python

from typing import List, Dict
from Bio import Entrez, SeqIO
import argparse
import urllib.request
from pathlib import Path
import pytaxonkit
import pandas as pd
from tqdm import tqdm
import sys

Entrez.email = "reni.denise@gmail.com"

######################################################
######################################################


def get_level_lineage(name: str) -> str:
    """Get level lineage from name .

    Args:
        name (str): The name of the level you want to know the lineage

    Returns:
        str: the level of the lineage
    """

    level_key = {
        "root": ["Viruses"],
        "realm": ["viria", "satellitia", "viroidia", "viriformia"],
        "subrealm": ["vira", "satellita", "viroida", "viriforma"],
        "kingdom": ["virae", "satellitae", "viroidae", "viriformae"],
        "subkingdom": ["virites", "satellitites", "viroidites", "viriformites"],
        "phylum": ["viricota", "satelliticota", "viroidicota", "viriformicota"],
        "subphylum": [
            "viricotina",
            "satelliticotina",
            "viroidicotina",
            "viriformicotina",
        ],
        "class": ["viricetes", "satelliticetes", "viroidicetes", "viriformicetes"],
        "subclass": [
            "viricetidae",
            "satelliticetidae",
            "viroidicetidae",
            "viriformicetidae",
        ],
        "order": ["virales", "satellitales", "viroidales", "viriformales"],
        "suborder": ["virineae", "satellitineae", "viroidineae", "viriformineae"],
        "family": ["viridae", "viriformidae", "viroidae", "satellitidae"],
        "subfamily": ["virinae", "satellitinae", "viroidinae", "viriforminae"],
        "genus": ["virus", "viriform", "viroid", "satellite"],
    }

    for level, exts in level_key.items():
        for ext in exts:
            if name.endswith(ext):
                return level
    return None


######################################################


def fix_taxa_column(lineage: List[str], species_name: str) -> str:

    """Return a string with prefix_taxa .

    Args:
        lineage (LIst(str)): List of all the lineage for a species
        species_name (str): Name of the species

    Returns:
        str: The full lineage with prefixes
    """

    prefix_taxa = {
        "root": "ro__",
        "realm": "r__",
        "subrealm": "sr__",
        "kingdom": "k__",
        "subkingdom": "sk__",
        "phylum": "p__",
        "subphylum": "sp__",
        "class": "c__",
        "subclass": "sc__",
        "order": "o__",
        "suborder": "so__",
        "family": "f__",
        "subfamily": "sf__",
        "genus": "g__",
        "species": "s__",
    }

    level = "root"  # If taxonomy empty there is a default value

    for name in lineage:
        if "unclassified" in name:
            level = False
        if "atellites" in name:
            level = False
        elif " viruses" in name:
            level = False
        elif " viroids" in name:
            level = False
        elif " phages" in name:
            level = False
        elif "unclassified Double-stranded satellite RNAs" in name:
            level = False
        elif "arenaviruses" in name:
            level = False
        elif (len(name.split()) > 1) or name == species_name:
            level = "species"
        else:
            level = get_level_lineage(name)
            try:
                level = level.rstrip()
            except:
                sys.exit(f"{lineage}")

        if level:
            prefix_taxa[level] += name

    if (level == "species" or level == "genus") and prefix_taxa["species"].endswith(
        "_"
    ):
        prefix_taxa["species"] = f"s__{species_name}"

    return ";".join(list(prefix_taxa.values()))


######################################################


def get_missing(all_not_there: List[str]) -> Dict[str, str]:

    """Get all missing entries for a genome.

    Args:
        all_not_there (List[str]): List of accession numbers

    Returns:
        Dict[str, str]: A dictionnary with all the missing informations and empty if not found
    """

    for acc_num in tqdm(
        all_not_there,
        colour="blue",
        desc="Missing genome done",
    ):

        missing = {
            "Representative": [],
            "Neighbor": [],
            "Host": [],
            "Taxonomy_name": [],
            "TaxID": [],
            "Lineage": [],
            "Lineage_prefix": [],
            "Rank": [],
            "Segment_name": [],
        }

        handle = Entrez.efetch(
            db="nucleotide", id=acc_num, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")

        taxa = ";".join(record.annotations["taxonomy"])
        missing["Lineage"].append(taxa)
        missing["Representative"].append(record.annotations["accessions"][0])
        missing["Neighbor"].append("")

        if record.features[0].qualifiers.get("host"):
            missing["Host"].append(record.features[0].qualifiers.get("host")[0])
        else:
            missing["Host"].append("")

        missing["Taxonomy_name"].append(record.annotations["organism"])

        if record.features[0].qualifiers.get("db_xref"):
            for elmt in record.features[0].qualifiers.get("db_xref"):
                if "taxon" in elmt:
                    taxon = elmt.split(":")[-1]
            missing["TaxID"].append(taxon)
        else:
            missing["TaxID"].append("")

        missing["Rank"].append("species")

        if record.features[0].qualifiers.get("segment"):
            segment = record.features[0].qualifiers.get("segment")[0]
            missing["Segment_name"].append(f"segment {segment}")
        else:
            missing["Segment_name"].append("segment")

        missing["Lineage_prefix"].append(
            fix_taxa_column(
                record.annotations["taxonomy"],
                species_name=missing["Taxonomy_name"][-1],
            )
        )
        missing["Lineage"][-1] = ";".join(
            name.split("_")[-1] for name in missing["Lineage_prefix"][-1].split(";")
        )

    return missing


######################################################


def get_metadata(refseq_tmp: Path, taxonomy_database: str) -> pd.DataFrame:

    """Read the refseq table into a pandas DataFrame

    Args:
        refseq_tmp (Path): Path to the refseq metadatafile downloaded
        taxonomy_database (str): Path to the taxonomy database folder

    Returns:
        pd.DataFrame: The dataframe with the addition of the lineage
    """

    refseq_df = pd.read_table(
        refseq_tmp,
        names=[
            "Representative",
            "Neighbor",
            "Host",
            "Selected lineage",
            "Taxonomy name",
            "Segment name",
        ],
        comment="#",
    )

    refseq_df["Name"] = refseq_df["Selected lineage"].str.split(",").str[-1]

    list_taxo_name = refseq_df["Name"].unique().tolist()

    # Change the name to taxonomic identifier
    name_taxid = pytaxonkit.name2taxid(list_taxo_name, data_dir=taxonomy_database)

    list_TaxID = name_taxid.TaxID.unique().tolist()

    # Change taxonomic identifier to lineage
    taxid_lineage = pytaxonkit.lineage(list_TaxID, data_dir=taxonomy_database)

    refseq_df = refseq_df.assign(
        Representative=refseq_df["Representative"].str.split(",")
    )
    refseq_df = refseq_df.explode("Representative")

    refseq_df = refseq_df.merge(name_taxid, on="Name", how="left")
    refseq_df = refseq_df.merge(
        taxid_lineage[["TaxID", "Lineage", "FullLineage"]], on="TaxID", how="left"
    )

    refseq_df = refseq_df.rename(
        columns={
            "Lineage": "Lineage_prefix",
            "FullLineage": "Lineage",
            "Selected lineage": "Selected_lineage",
            "Taxonomy name": "Taxonomy_name",
            "Segment name": "Segment_name",
        }
    )

    refseq_df = refseq_df[
        [
            "Representative",
            "Neighbor",
            "Host",
            "Taxonomy_name",
            "TaxID",
            "Lineage",
            "Lineage_prefix",
            "Rank",
            "Segment_name",
        ]
    ]

    refseq_df = refseq_df.drop_duplicates().reset_index(drop=True)

    return refseq_df


######################################################


def parse_args():

    """parse command line arguments

    Returns:
        argparse.parser : All the argument parsed
    """
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This script will create the metadata file with taxonomy annotation for each genome in RefSeq"
    )

    # Positional mandatory arguments
    parser.add_argument(
        "--refseq_nucl_file", help="The refseq nucleotide file", type=str, required=True
    )

    parser.add_argument(
        "--refseq_version",
        help="The refseq version (e.g 211, 216...)",
        type=int,
        required=True,
    )

    parser.add_argument(
        "--output_folder",
        help="The name of the output folder, default: the same folder as the fasta file",
        type=str,
        default="",
    )

    parser.add_argument(
        "--taxonomy_database",
        help="Path to the NCBI taxonomy database folder",
        type=str,
        required=True,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


######################################################
######################################################


def main(args):

    OUTPUT_FOLDER = (
        Path(args.output_folder)
        if args.output_folder
        else Path(args.refseq_nucl_file).parent
    )
    refseq_tmp = OUTPUT_FOLDER / f"refseq_{args.refseq_version}_metadata.tsv"
    refseq_tmp_torm = OUTPUT_FOLDER / f"refseq_{args.refseq_version}_metadata.tmp.tsv"
    refseq_output = OUTPUT_FOLDER / f"refseq_{args.refseq_version}_metadata.lineage.tsv"

    url = "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"
    file_name = refseq_tmp

    message = "Downloading the metadata from NCBI"

    print("-" * (len(message) + 2))
    print(f"| {message}")
    print("-" * (len(message) + 2))

    if not file_name.exists():
        response = urllib.request.urlretrieve(url, file_name)

    print("\n Done!\n")

    message = "Getting the taxonomy lineage"
    print("-" * (len(message) + 2))
    print(f"| {message}")
    print("-" * (len(message) + 2))

    if refseq_tmp_torm.exists():
        refseq_df = pd.read_table(refseq_tmp_torm)
    elif not refseq_output.exists():
        refseq_df = get_metadata(
            refseq_tmp=refseq_tmp, taxonomy_database=args.taxonomy_database
        )

        records = refseq_df.to_dict("records")
        good_lineage = []

        for row in tqdm(
            records,
            colour="blue",
            desc="Adding prefix to lineage",
        ):
            good_lineage.append(
                fix_taxa_column(
                    lineage=row["Lineage"].split(";"), species_name=row["Taxonomy_name"]
                )
            )

        # Adding the good lineage
        refseq_df["Lineage_prefix"] = good_lineage

        # Make sure that all the lineage is there
        refseq_df["Lineage"] = refseq_df["Lineage_prefix"].apply(
            lambda row: ";".join([name.split("_")[-1] for name in row.split(";")])
        )

        refseq_df.to_csv(refseq_tmp_torm, index=False, sep="\t")

    print("\n Done!\n")

    message = "Getting missing information"
    print("-" * (len(message) + 2))
    print(f"| {message}")
    print("-" * (len(message) + 2))

    if not refseq_output.exists():
        list_taxo_ncbi = set(refseq_df.Representative.tolist())

        list_taxo_nucl = set(
            [
                genome.id.split(".")[0]
                for genome in SeqIO.parse(args.refseq_nucl_file, "fasta")
            ]
        )

        all_not_there = list_taxo_nucl - list_taxo_ncbi

        missing = get_missing(all_not_there=all_not_there)

        refseq_df = pd.concat([refseq_df, pd.DataFrame(missing)], ignore_index=True)

        refseq_df.to_csv(refseq_output, index=False, sep="\t")

        refseq_tmp_torm.unlink()

    print("\n Done!\n")

    return


######################################################
######################################################


if __name__ == "__main__":
    args = parse_args()
    main(args)
