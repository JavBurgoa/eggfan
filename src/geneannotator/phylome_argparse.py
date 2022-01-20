from typing import Annotated
import phylome
import utils
import pandas as pd
import os
import argparse
import shutil


# 1. which method?
#     A. HGNC?
#     annotated_tables = phylome.annotate_orthology_HGNC_method(args.query, args.ortho_tables)

#     B1. read or make lookup?
#     lookup = phylome.make_lookup(args.ortho_tables)

#     B2. save lookup?

#     B3. read or make translated_orthologies?
#     translated_orthologies = phylome.translate_orthologies(args.ortho_tables, lookup)

#     B4. save translated_orthologies?

#     B5. run rest
#     annotated_tables = phylome.find_query_orthologs(args.query, translated_orthologies)

# 2. save result
# phylome.save_annotated(annotated_tables, args.output, args.suffix)

# flags = [HGNC]


def main(query, ortho_tables, output, input_lookup, suffix, flags):

    if suffix == None:
        suffix = "_annotated_orthology"

    # first go to the result directory
    #os.chdir(output)
    if flags["HGNC"]:
        annotated_tables = phylome.annotate_orthology_HGNC_method(query, ortho_tables)
    else:
        lookup = get_lookup(ortho_tables, input_lookup)
        translated_orthologies = get_translated_orthologies(
            ortho_tables, lookup, flags["input_translated"]
        )
        annotated_tables = phylome.find_query_orthologs(query, translated_orthologies)
    
    print("saving") # These two lines below save as long as you didn't input the lookup and/or the translated tables
    save_lookup(lookup, output, flags["input_translated"], input_lookup)
    save_translated(translated_orthologies, output, flags["input_translated"])
    
    phylome.save_annotated(annotated_tables, output, suffix)
    print("done")








def get_translated_orthologies(ortho_tables, lookup, input_translated):
    """
    either read the orthology table(s) or make them
    """
    if input_translated:

        translated_orthologies = phylome.read_translated_tables(ortho_tables)

    else:
        
        translated_orthologies = phylome.translate_orthologies(
                ortho_tables, lookup
            )
    return translated_orthologies


def get_lookup(ortho_tables, lookup, overwrite=False):
    if lookup is not None:
        lookup = pd.read_csv(lookup, sep="\t", keep_default_na=False)
    else:
        lookup = phylome.make_lookup(ortho_tables)
    return lookup


def read_translated_ontologies(path_to_translated_orth):
    translated_orthologies = list(pd.read_csv(path_to_translated_orth, sep="\t"))
    return translated_orthologies


def save_lookup(lookup, output, input_translated, input_lookup):
    """
    Save lookup or not depending on context
    """
    no_save = [not input_translated, not isinstance(input_lookup, str)]
    if all(no_save): # Only if there is no input_trans and no input_lookup then save
        lookup.to_csv(output + "lookup.tsv", sep="\t")


def save_translated(annotated_tables, output, input_translated):
    """
    Save translated orthotables or not depending on context
    """
    if not input_translated:
        save_dir = output + "/translated_orthology_tables"
        if os.path.exists(save_dir):
            shutil.rmtree(save_dir)
        os.mkdir(save_dir)
        phylome.save_annotated(annotated_tables, save_dir, suffix = "_translated")




if True:
    ##### Arguments parser #####
    parser = argparse.ArgumentParser(
        description="Subsets phylome orthology tables to include only proteins that have as orthologs genes/proteins from your human query. Then adds somecolumns to excplicit which genes/proteins are the orthologs"
    )

    # Full pipeline
    parser.add_argument(
        "-t",
        "--ortho-tables",
        type=str,
        dest="ortho_tables",
        metavar="DIR",
        help="Path to folder containing all the unadulterated phylome orthology files you want to annotate or a path to one of those files.",
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        metavar="FILE",
        help="Path to file with curated list of human genes representing you module/family of interest. Genes should be in ENSEMBL_GENE_ID format with default options or HGNC (genecards) format with --HGNC flag",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="",
        required=True,
        help="Path to folder where you want to save outputs. Can be saving a lookup, translated tables or annotated tables.",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        metavar="",
        required=False,
        help="Optional. name of output files will be <taxID><suffix>.tsv . Default '_annotated_orthology' for annotated files, '_translated' for translated files. Use it if you are outputting the annotated files (end of pipeline) or --output_translated_orthotables",
    )

    parser.add_argument(
        "-l",
        "--lookup",
        type=str,
        metavar="",
        required=False,
        help="Optional. Path to lookup table that will be used to make translated orthology tables"
    )

    parser.add_argument(
        "--input_translated",
        action="store_true",
        help="If true a new translation of the orthology tables will not be made. The pipeline will be run with the orthology tables in --ortho_tables ",
    )
    # HGNC method
    parser.add_argument(
        "--HGNC",
        action="store_true",
        dest="hgnc",
        help="Use HGNC to do the matching. This avoids making a lookup and translating the orthology tables",
    )

    args = parser.parse_args()
    flags = {}
    flags["input_translated"] = args.input_translated
    flags["HGNC"] = args.hgnc

    if __name__ == '__main__':
       main(args.query, args.ortho_tables, args.output, args.lookup, args.suffix, flags)
