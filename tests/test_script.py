from eggfan import phylome
from eggfan import goterms
from eggfan import orthogroup
import pandas as pd

# GOTerms
cap_annotated = goterms.GOTerms_annotation("tests/data/Capitella_emapper_redux.txt", GOterm="GO:0003700")
print(cap_annotated)
exit()


# Orthogroup

eggnog = orthogroup.read_eggnog(['tests/data/eggnog/Eggnog_Bilateria(33213)_members.tsv', 'tests/data/eggnog/Eggnog_Metazoa(33208)_members.tsv'])
lookup = pd.read_csv('tests/data/eggnog/lookup.txt', sep='\t')
query = pd.read_csv("tests/data/TFs_Human_Ensembl.txt", sep = "\t")
emapper = pd.read_csv("tests/data/Capitella_emapper_redux.txt", skiprows=4, sep="\t")

translated_eggnog = orthogroup.egg_translate(eggnog, lookup)
orthogroup.translated_QC(translated_eggnog)
orthogroup.translated_QC(translated_eggnog, query)
query_orthogroups = orthogroup.merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
annotated_genes = orthogroup.emapper_annotation(emapper, query_orthogroups, keep_all_targets= False)




# Phylome regular version

lookup = phylome.make_lookup("tests/7227_orthologs_redux_1000.tsv")
print(lookup)
lookup = pd.read_csv("tests/data/lookup.tsv", sep = "\t")
print(lookup)

translated_orthologies = phylome.translate_orthologies("tests/data/phylomes/7955_Zebrafish_orthologos_redux.tsv", lookup)

print(translated_orthologies)

phylome.save_annotated(translated_orthologies, "tests/data/translated_orthotables/", suffix = "_translated")

annotated_tables = phylome.find_query_orthologs("tests/data/TFs_Human_Ensembl.txt", translated_orthologies)

phylome.save_annotated(annotated_tables, "tests/results/")


## Phylome HGNC version


annotated_tables = phylome.annotate_orthology_HGNC_method("tests/data/TF_Human_HGNC.txt", "tests/data/translated_orthotables/")
phylome.save_annotated(annotated_tables, "tests/results/")

