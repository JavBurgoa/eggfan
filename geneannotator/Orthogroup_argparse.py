import orthogroup
import argparse
import pandas as pd
import numpy as np

##### Arguments parser #####
parser = argparse.ArgumentParser(description="Gives you all genes/proteins from emapper annotated based on a list of curated human genes")
parser.add_argument("-eg", "--eggnog", action='append', type=str, metavar="",required=True, help="one or several paths to eggnog members files (de-compressed). For each path add an '-eg' flag")
#use like: function --eggnog "path" --eggnog "path" -eg "path"
parser.add_argument("-l", "--lookup", type=str, metavar="",required=True, help="Path to lookup table for eggnog translation")
parser.add_argument("-q", "--query", type=str, metavar="",required=True, help="Path to file with curated list of genes")
parser.add_argument("-em", "--emapper", type=str, metavar="",required=True, help="Path to output file from emapper, with all target proteins to annotate")
parser.add_argument("-m", "--merge_on", type=str, metavar="",required=True, help="Name of column from lookup file that will be ued to make the matching. This is, the name of the column with the final translation of the genes")
parser.add_argument("--QC", action="store_true", help="Print QC of eggnog translation and stop pipeline")
parser.add_argument("--keep_all_targets", action="store_true", help="keep all genes from emapper whether they match with the query list or not")
parser.add_argument("--remove_conversions", action="store_true", help="Remove all gene ID conversion from the lookup table, keep only the IDs used in matched_column")
args = parser.parse_args()


#### Pipeline #####
## Load datasets
eggnog = orthogroup.read_eggnog(args.eggnog)
lookup = pd.read_csv(args.lookup, sep='\t')
query = pd.read_csv(args.query, sep = "\t")
emapper = pd.read_csv(args.emapper, skiprows=4, sep="\t")


## Run pipeline
# translate eggnog Protein ENSEMBL IDs to whatever you want (default and recommended, ENSEMBL gene IDs)
translated_eggnog = orthogroup.egg_translate(eggnog, lookup)
#print("* Eggnog files read")

# QC
if args.QC:
	orthogroup.translated_QC(translated_eggnog, query)
	exit("QC finished, if you want to run the full pipeline remove the flag '--lookup'")

# Make table of all your query genes and their respective eggnog orthogroups.
if args.remove_conversions:
	keep_conversions = False
else: keep_conversions = True

query_orthogroups = orthogroup.merge_with_query(translated_eggnog, query, merge_on = args.merge_on, keep_conversions= keep_conversions)
#print("* All files read and Query - Orthogroup table made")

# Match my query genes' orthorgoups with my proteome emapper orthogroups
if args.keep_all_targets:
	keep_all_targets = True
else: keep_all_targets = False
annotated_genes = orthogroup.emapper_annotation(emapper, query_orthogroups, keep_all_targets)


# Output
print(annotated_genes.to_csv(sep='\t', index=False))