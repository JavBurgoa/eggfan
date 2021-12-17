import phylome
import pandas as pd
import os

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


def main(query, ortho_tables, output, suffix, flags):
    # first go to the result directory
    os.chdir(output)
    if flags.HGNC:
        annotated_tables = phylome.annotate_orthology_HGNC_method(query, ortho_tables)
    else:
        lookup = get_lookup(ortho_tables)
        translated_orthologies = get_translated_orthologies(
            ortho_tables, lookup, flags.overwrite
        )
        annotated_tables = phylome.find_query_orthologs(query, translated_orthologies)
    phylome.save_annotated(annotated_tables, output, suffix)


def get_translated_orthologies(ortho_tables, lookup, overwrite=False):
    species_id = phylome.get_species_id(ortho_tables)
    translated_path = "./translated/" + species_id + "_human_orthologs.tsv"
    if os.path.isfile(translated_path):
        translated_orthologies = list(pd.read_csv(translated_path, sep="\t"))
    else:
        translated_orthologies = phylome.translate_orthologies(
            ortho_tables, lookup, overwrite=overwrite
        )
    return translated_orthologies


def get_lookup(ortho_tables, overwrite=False):
    if os.path.isfile("./lookup.tsv"):
        lookup = pd.read_csv("./lookup.tsv", sep="\t", keep_default_na=False)
    else:
        lookup = phylome.make_lookup(ortho_tables, overwrite=overwrite)
    return lookup


def read_translated_ontologies(path_to_translated_orth):
    translated_orthologies = list(pd.read_csv(path_to_translated_orth, sep="\t"))
    return translated_orthologies
