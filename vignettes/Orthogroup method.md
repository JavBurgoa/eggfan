# **Orthogroup Method**
This method uses the following basis: If you know a gene has a particular function in humans (The gene is a transcription factor, or a contractile module, etc.) and another gene (target) from another species belongs to the same orthogroup than that human gene, then the target gene also probably is a transcription factor, contractile module etc in its species.

Below you will find a short tutorial on how to use this pipeline in python. We will use the following example: We want to know which genes in the *Capitella* *teleta* genome are transcription factors (TFs).

&nbsp;
&nbsp;

## **Required files**
There are some files we will need in order to run this pipeline. The current version 0.0.1 (Oct-2021) needs to gather some files manually, future versons will automatize this procedure.

### Query list
To run the pipeline, we need a list with human genes that represent all genes with the desired function in humans, in this case, a curated list with all human TFs. This information is present in published articles, if you want to test a different gene module you'll have to find the information by yourself. Once you have gathered the information simply wright a text/tsv/csv file with a single line per gene with no header. Mandatory Ensembl Gen ID format. You can find a sample in test/data/TFs_Human_Ensembl.txt

&nbsp;
&nbsp;
### Eggnog database
We need to know for each of the genes in the "query" list, what are their respective orthogroups. This data is stored in the Eggnog database. We can use orthogroups at different levels. Let's say we want to relate Human and *Capitella* genes by both Metazoan and Bilaterian orthogroups. Then we need to go to the [Eggnog database](http://eggnog5.embl.de/#/app/home), go to Downloads, find out our level(s) of interest, in this case Metazoan and Bilateria, and download the (taxID)\_members.tsv.gz file one at a time. De-compress it(them) and save in a folder.

&nbsp;
&nbsp;

### Lookup for conversions
The eggnog files we just downloaded show us all orthogroups witin a tax level and all **Proteins** that belong to that group. However, in most cases our query won't be a list of **Proteins** but a list of **Ensembl Gen IDs** or something else. So we need to translate those Proteins in the eggnog files to whatever we have in our query list. To do so we need a lookup table that gives us for each protein in eggnog, its Ensembl Gen ID equivalent.

For this example you can use the table in test/data/lookup.tsv, but you should ownload your own from Biomart for things other than TFs.

The recommended lookup table would be human Ensembl protein IDs to Ensembl gene IDs and HGNC symbols (this last one is optional but highly recommended to add). You can find an example of this lookup in /test/data/lookup.tsv in the github repository.

> However you could have your query list in something other than Ensemble gene IDs, in which case you could switch the Enseble gene IDs column of the lookup table by whatever ID
> format you have in your query. However this hasn't been very well tested and may lead to problems. At the moment this pipeline only works if your query contains human genes.

To do so we go to [Biomart](https://www.ensembl.org/biomart/martview/62f7cff64a6aaf9711bf7c0a3e52f7e7), and we select the species we want to translate, in this case, human. Then, in attributes we choose what are the conversions we want to make. We are hoosing **Protein stable** ID and **Gene stable ID**. Additionally, we want to know a more human readable name for each gene, so we also go to **EXTERNAL** and select HGNC symbol. Finally we click on Result in the top left of the screen and after it loads we hit GO! to download in tsv format.

> You could run the pipeline without the translation, for which, at the moment (Oct-2021)
> there is no in-built functionality, however, you can fake this by creating
> a lookup table where there are two identical columns, both with Protein IDs

&nbsp;
&nbsp;

### Emapper output
Finally we need to know to which orthogroups each of our *Capitella* target genes belong to, in order to match this with the Eggnog database information. 
To do so, we need a proteome (or genome) in FASTA format, in our case *Capitella teleta*'s proteome. We run it using [emapper](http://eggnog-mapper.embl.de/) with standard values. We save the output file in a folder as is.

&nbsp;
&nbsp;

## **The pipeline**
### Install package
To do so follow the README.md file in our GitLab repository

### **Python**
### Import modules
On the top of your file put the required dependencies:
```
>>> from eggfan import orthogroup
>>> import pandas as pd
```

#### Import data
First we need to load all data we mentioned in the [Required data](##required-files) section. We will show now how to import all these datasets.

1. **Eggnog datasets**
For this we have a function called read_eggnog, that allows you to import as many tax levels as you decide. Simply add as arguments a path for every (de-compressed) eggnog file you have. Remember that if you have more than one eggnog file the path should come inside of a list [].
```
>>> eggnog = read_eggnog(['tests/data/eggnog/Eggnog_Bilateria(33213)_members.tsv', 'tests/data/eggnog/Eggnog_Metazoa(33208)_members.tsv'])
```

2. **Lookup table**
Simply import it using the pandas function read_csv. Don't forget to add whatever separator the file has in the *sep* argument ("\t" for tsv files, "," for csv, etc.). For this example you can use the table in test/data/lookup.tsv
```
>>> lookup = pandas.read_csv('/test/data/lookup.tsv', sep='\t')
```

3. **List of Human transcription factors**
Again, just a text/tsv/csv file with a single column containing all your curated genes. In our case, human transcription factors. Put a header in the file (some text on the first line, a name for the column) or specify "header = None". For TFs we left a list with all human TFs in test/ folder.
```
>>> query = pd.read_csv("test/data/TFs_Human_Ensembl.txt", header=None, sep = "\t")
```

4. **Emapper output**
Target genes emaper output exactly as it comes from emapper. The argument *skiprows=4* removes the top 4 row of the file, containing metadata information. Do not remove this argument. For Capitella we left a emapper output in test/ folder.
```
>>> emapper = pd.read_csv("test/data/Capitella_emapper_redux.txt", skiprows=4, sep="\t")
```

&nbsp;
&nbsp;

### Run pipeline

First we need to create a table with all our human transcription factors and their respective orthogroups at all the taxonomic levels chosen. This can be breaken in two parts: Translating the eggnog database from protein IDs to gene IDs and then removing from there all genes that are not in our query (list of human TFs).

First we translate the eggnog datasets with egg_translate() like this
```
>>> translated_eggnog = orthogroup.egg_translate(eggnog, lookup)
```

Then we can run a quality control script to tell us wether any of the genes in the eggnog database were not translated
```
>>> orthogroup.translated_QC(translated_eggnog)
```

If there were any non-translated genes, then the lookup table didnot catch all eggnog proteins, wether because they are outdated or anything else. 
We can check if any of the non-translated genes were not translated. If you find any of your curated genes is outputted by this function
you should go check in the eggnog database wether this isnbecause your gene wasn't in the database at all (in which case, nothing should be done) or your gene was in the eggnog database but was not translated, in which case the lookup should be updated to account for that otherwise lost gene.
```
>>> orthogroup.translated_QC(translated_eggnog, query)
```

Back to the pipeline, once we have checked the translation was optimal we can proceed by removing from the translated output any gene that is not a transcription factor. We use the following function. In merge_on we put the name of the column of translated_eggnog where we want to do the merge.
```
>>> query_orthogroups = orthogroup.merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
```

This outputs a table with all your TFs, and their respective orthogroups at the levels chosen like this
|Gene       |Orthogroup@33208|Orthogroup@33298|HGNC        |ProteinID      |
| ------ | ------ | ------ | ------ | ------ |
|Query Gene1|FGFV@33208      |BVGC@33298      |Gene symbol1|Target Protein1|
|Query Gene2|TGOK@33208      |CACF@33298      |Gene symbol2|Target Protein2|

> The keep_conversions argument allows the user to remove any other gene conversions 
> introduced by the lookup table. In our case it would remove the HGNC and ProteinID columns

The Orthogroup@ columns show you the orthogroups at the different taxonomical levels. The level is denoted by the @(taxID)

If we combine that with the emapper output we can know which of the *Capitella* genes (emapper) share an orthogroup at Metazoan or Bilaterian level with a Human TF. We can do this with the following code:
```
>>> annotated_genes = orthogroup.emapper_annotation(emapper, query_orthogroups, keep_all_targets= False)
```

annotated_genes contains a table with all *Capitella* genes that are putative TFs, based on your list of Human TFs and the taxonomic levels chosen at Eggnog. You can save this file as a tab separated file with the following code:
```
>>> annotated_genes.to_csv("/g/arendt/Javier/Python/eggfan/tests/Ortho_method_Capitella_TFs.tsv", sep = "\t")
```

&nbsp;
&nbsp;

### **Command line**
If you want to run the orthogroup pipeline directly from command line you don't need to import the data as specified before. You simply need the paths to your files.
However, you will need to activate python in your commandline. You can check if you have python already available running `python -v`. If you don't have it installed a vey easy way is with a Miniconda environment. However, this is not explaned in this vignette.

You can see all options in the command line orthogroup function simply with
```
python src/eggfan/orthogroup_argparse.py -h
```

#### Run pipeline
There are several mandatory flags needed to run this program:
1. **-g** 
You need to put here a path or several paths to the eggnog dasets as specified in [Eggnog database](###eggnog-database). If you want to put more than one path separate each path with a new flag
2. **-l**
Path to a lookup table as described in [Lookup for conversions](###lookup-for-conversions).
3. **-q**
Path to file with a single column containing all your human query genes [Query list](###-query-list).
4. **-e**
Path to emapper output run with your target proteins you want to annotate [Emapper output](###-emapper-output).
5. **-m**
Name of the column from the lookup table with the ID format in which your query is.

With this information you can now run the pipeline from command line. Here is an example usage for the abe example as with thepython version: Which genes in the *Capitella* *teleta* genome are transcription factors (TFs).

```
python eggfan/Orthogroup_argparse.py -g 'tests/data/eggnog/Eggnog_Bilateria(33213)_members.tsv' -g 'tests/data/eggnog/Eggnog_Metazoa(33208)_members.tsv' -l 'test/data/lookup.txt' -e "test/data/Capitella_emapper_redux.txt" -q "test/data/TF_Human_Ensembl.txt" -m "Gene stable ID"
```
This will output a dataframe with all of your emapper genes that share orthogroup with your human query genes and their respective translations. the same functionality than the python version, but in a single line of code.

There are also some additional flags one can add to modify the output. For example `--keep_all_targets` allows you to have all genes in the emapper file in the inal output, whether they share orthorgoup with your query or not.

`--rm_conversions` leaves you only with the list of annotated genes, with no translated versions of them.

Finally `--QC` does not modify the output. Actually it doesn't even let the pipeline finish. It outputs a list with all proteins in eggnog that were not translated. This can happen because the lookup table is defective or because some of eggnog's protein IDs are outdated and therefore do not match the lookup table. You can correct this by adding the outdated versions manually in the lookup.