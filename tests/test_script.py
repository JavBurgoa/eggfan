from geneannotator import phylome
import pandas as pd

lookup = phylome.make_lookup("/g/arendt/Javier/Python/test_geneannot/geneannotator/tests/7227_orthologs_redux_1000.tsv")
print(lookup)
lookup = pd.read_csv("/g/arendt/Javier/Python/test_geneannot/geneannotator/tests/test1/lookup.tsv", sep = "\t")
print(lookup)

translated_orthologies = phylome.translate_orthologies("/g/arendt/data/phylomeV2/orthology_tables_nocollapse/7955_orthologs.tsv", lookup)

print(translated_orthologies)

phylome.save_annotated(translated_orthologies, "tests/test3", suffix = "_translated")

annotated_tables = phylome.find_query_orthologs("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", translated_orthologies)

phylome.save_annotated(annotated_tables, "tests/test2/")