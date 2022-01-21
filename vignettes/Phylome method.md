# **Phylome method**
This method does not rely on Eggnog or emapper directly. Instead it works on some of the output files from the Phylome, another resource created by the Jaime Huerta Cepas lab in CBGP, Madrid, Spain.
The phylome uses proteomes from a great number of species. One of the outputs from the phylome are orthology tables that show for each of the genes in each of those species, which genes are their orthologs.
This tables show, for each species' gene, what are their orthologs (in all the species in the phylome). As humans are in the phylome, these tables tell you all huaman orthologs of any species in the phylome.

Having a list of human genes representing a gene module/family we can know for any gene in the phylome, if they are orthologous to a human family/module. If they are we can say the target gene also belongs to the same module/family (putatively).

For demonstration purposes we will now run a pipeline where we annotate genes in two random phylome species with the question:  Which genes in the *Cladocora caespitosa* and *Danio rerio* (Zebrafish) proteomes are transcription factors (TFs)?.

## **Required files**
### Phylome orthology tables
Non-collapsed orthology tables outputted from the phylome. They should be in a foder **containing only orthology tables**. If you just want to annotate few species you can keep a copy of those in a separate folder. They will not be overritten.

### Query list
To run the pipeline, we need a list with human genes that represent all genes with the desired function in humans, in this case, a curated list with all human TFs. To do so, simply wright a text/tsv/csv file with a single line per gene. GeneIDs can be in either EnsemblID shape for the regular pipeline or HGNC (Genecards) format or the HGNC pipeline.


## **The pipeline**
There are two versions of the pipeline, the regular and the HGNC.
The orthlogy tables output orthologs in HGNC and Uniprot format. If you want to match your query list with the orthology tables using HGNCs (**HGNC method**), then there is only one step for the pipeline. If you want to do the matching with ENSEMBLIDs (**regular method**) you need to first create a lookup table and then use it to translate the orthology tables to EnsemblIDs and finally fo the matching, 3 steps.

### Python
### Install package
To do so follow the README.md file in our GitLab repository

### Import modules
For this pipeline you don't need to import any file. Just open your script, and import the following packages
```
>>> from geneannotator import phylome
>>> import pandas as pd
```

### Import data
You don't need to import any data if it's the first time you are running the pipeline or you are running the HGNC method.
If however you have already run the pipeline before and you have a lookup table you can just import it with:

```
>>> lookup = pd.read_csv("path/to/lookup.tsv", sep = "\t", keep_default_na=False)
```

If you are running the regular pipeline and you already have the orthology tables translated to ENSEMBLID we will see later how to use them.


### Run pipeline

1. **Make a lookup**
If it's the first time running the pipeline first we have to make a lookup table to translate the orthology table(s) to ENSEMBLIDs. This procedure makes a lookup with the human genes that appear in the inputted tables. The recommended behaviour is saving the lookup table and the subsequent trandlations of the phylome. This way you avoid making them every time as it is quite time consuming.
```
>>> lookup = phylome.make_lookup("path/to/folder/with/orthology_tables")
```
> Alternatively you can make a lookup for a few or a single orthology table by specifying having only few species in the folder or putting a path to a single file. 
> However, this is  not recommended

2. **Save the lookup**
Now you can save the lookup table like:
```
>>> lookup.to_csv("save/path/lookup.txt, sep = "\t", index = False)
```
Alternatively you can not save it and continue with the pipeline

3. **Translate orthology table(s)**
Now you need to translate the orthology tables human orthologs from UniprotKB to ENSEMBL gene ID. To do so use the newly created lookup (or the one you saved, see [Import data](###import-data)) 
```
>>> translated_orthologies = phylome.translate_orthologies("path/to/phylome/orthology/tables", lookup)
```
You can also just input a path to a single orthology table.
If you want to save these to be able to do multiuple annotations on them without running the translation again you can do so with the following code:

```
>>> phylome.save_annotated(translated_orthologies, "output/directory/", suffix = "_translated")
```
The last argument is the suffix. the name of output files will be taxID_suffix.tsv. You can change the suffix by changing that parameter.

4. **Annotate genes**
Finally to annotate genes simply use the translated tables and your human query (gene family/module).
```
>>> annotated_tables = phylome.find_query_orthologs("/path/to/query", translated_orthologies)
```

If you already had the translated tables you can just specify path to directory or file:
```
>>> annotated_tables = phylome.find_query_orthologs("/path/to/query", "path/to/translated_orthologies")
```

5. **Save annotated tables**
annotate_files is a list containing all your annotated datasets. To save them (recommended) use the following function:

```
>>> phylome.save_annotated(annotated_tables, "/path/to/saving/folder/")
```
Name of the output files will be taxID_annotated_orthology.tsv. You can change the "annotated_orthology" suffix with the `suffix` parameter.

The final tables will show all genes in the original orthology table that have as ortholog any human gene in your query like this:

|##Seed_(co-)orthologs|type|ENSEMBL_ID|orthologs|GeneName_target|Ensembl_query-only|GeneName_target_query-only|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
|Phylome Gene1|orthology type|All Human orthologs ENS_ID|All Human orthologs UNIPROT_ID|All Human orthologs HGNC|Query Human orthologs ENS_ID| Query Human orthologs HGNC|

- ##Seed_(co-)orthologs: Genes directly from phylome orthology tables. Only genes that have as ortholog a human gene from the query list
- type: orthology type directly from phylome orthology tables
- ENSEMBL_ID: ENSEMBLID conversion of ALL the human orthologs of ##Seed_(co-)ortholog
- orthologs: ALL the human orthologs of ##Seed_(co-)ortholog in Uniprot format. This column comes directly from phylome orthology tables.
- GeneName_target: HGNC conversion of ALL the human orthologs of ##Seed_(co-)ortholog. This column comes directly from phylome orthology tables.
- Ensembl_query-only: Gene(s) from ENSEMBL_ID column that are also in the query. The genes that are not in the query are represented with dashes.
- GeneName_target_query-only: Gene(s) from ENSEMBL_ID column that are also in the query translated to HGNC. The ones that are not are represented with dashes.


### HGNC version
This ipleine is simpler, faster and yields simmilar results than the regular pipeline. However the slight differences might be worth checking depending on the case.
There is no need to translate or create a lookup table. Simply, using your query and the phylome you can run:
```
>>> annotated_tables = phylome.annotate_orthology_HGNC_method("path/to/query.txt", "path/to/orthology_tables/")
```
Alternatively you can put a path to a single file.

Finally, save your annotated orthology tables with:
```
>>> phylome.save_annotated(annotated_tables, "/path/to/saving/folder/")
```
The output is the same as with the regular method but without the ENSEMBL_ID and Ensembl_query-only columns.


############################ Not like this anymore beloww VVVVV

### Command line
To run the pipeline in command line you have the `python src/geneannotator/phylome_argparse.py` function. It works exactly the same as the pure python version, but in a single line.

If you want to know what each flag does, just run on your terminal:
```
python src/geneannotator/phylome_argparse.py -h
```

You can run this single line:
```
python src/geneannotator/phylome_argparse.py -t "path/to/phylome_tables/" -q "path/to/human_query.tsv" -o "saving/path/"
```
> Remember using `--suffix` if you want to change the name of the final output(s) (the final names will be taxID_suffix.tsv).
This line will go through the following process:

1. **Save a lookup**
The pipeline will create a lookup table for all the human orthologs in the phylome files. This will be saved in the output directory.

2. **Translate**
Then it will translate the phylome orthology tables using our saved lookup. Then it will create a folder called "translated_orthology_tables" in the output directory. All the phylome files with their human orthogroups translated to ENSEMBLID in an extra column will be saved here. 

3. **Annotate**
Finally, it goes over the translated tables and annotates the genes based on orthology. This ouput is saved as separate files in the output directory.


### Re-running with a different query
You may want to annotate the same species with a different module. In that case you don't need to run the whole pipeline again. You can specify that you already have the translated tables with the `--input_translated`:
```
python src/geneannotator/phylome_argparse.py -t "path/to/translated_orthology_tables/" -q "path/to/human_query.tsv" -o "saving/path/" --input_translated
```
This way you bypass making the lookup and translations aswell as saving them.


### **HGNC method**
Much simpler than the regular method. It will use the already present in the orthology tables HGNCs to make the matchings. You only need your orthology table(s) and your query (in HGNC format), and then run:
```
python src/geneannotator/phylome_argparse.py -t "path/to/phylome/orthology_tables/folder/" -q "path/to/human_query.tsv" -o "path/to/output/folder" --HGNC_method
```
The output are your annotated orthology tables but this time without the ENSEMBL translation columns.
Remember using `--suffix` if you want to change the names of the outputs (the final names will be taxID_suffix.tsv).
