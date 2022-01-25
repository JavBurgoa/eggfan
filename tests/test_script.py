from geneannotator import phylome
from geneannotator import goterms
from geneannotator import orthogroup
import pandas as pd

# Phylome
"""
lookup = phylome.make_lookup("/g/arendt/Javier/Python/test_geneannot/geneannotator/tests/7227_orthologs_redux_1000.tsv")
print(lookup)
lookup = pd.read_csv("/g/arendt/Javier/Python/test_geneannot/geneannotator/tests/test1/lookup.tsv", sep = "\t")
print(lookup)

translated_orthologies = phylome.translate_orthologies("/g/arendt/data/phylomeV2/orthology_tables_nocollapse/7955_orthologs.tsv", lookup)

print(translated_orthologies)

phylome.save_annotated(translated_orthologies, "tests/test3", suffix = "_translated")

annotated_tables = phylome.find_query_orthologs("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", translated_orthologies)

phylome.save_annotated(annotated_tables, "tests/test2/")
"""
## HGNC version
"""
annotated_tables = phylome.annotate_orthology_HGNC_method("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TF_names_v_1.01.txt", "/g/arendt/Javier/Python/geneannotator/tests/translated_orthotables/")
phylome.save_annotated(annotated_tables, "/g/arendt/Javier/Python/geneannotator/tests/")
"""


# GOTerms
"""
cap_annotated = goterms.GOTerms_annotation("/g/arendt/Javier/Python/TF_annot_methods/Capitella_teleta/Data/Capitella_teleta_Gene_emapper_annotations.txt", GOterm="GO:0003700")
print(cap_annotated)
"""


# Orthogroup
"""
eggnog = orthogroup.read_eggnog(['/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Bilateria(33213)_members.tsv', '/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Metazoa(33208)_members.tsv'])
lookup = pd.read_csv('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')
query = pd.read_csv("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", sep = "\t")
emapper = pd.read_csv("/g/arendt/Javier/Python/TF_annot_methods/Capitella_teleta/Data/Capitella_teleta_Prot_emapper_annotations.txt", skiprows=4, sep="\t")

translated_eggnog = orthogroup.egg_translate(eggnog, lookup)
orthogroup.translated_QC(translated_eggnog)
orthogroup.translated_QC(translated_eggnog, query)
query_orthogroups = orthogroup.merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
annotated_genes = orthogroup.emapper_annotation(emapper, query_orthogroups, keep_all_targets= False)
"""