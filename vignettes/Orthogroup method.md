# Orthogroup Method
Ths method uses the following basis: If you know a gene exerts a particular function in humans (The gene is a transcription factor, or a contractile module, etc.) and another gene (target) from another species belongs to the same orthogroup than that human gene, then the target gene also probably is a transcription factor, contractile module etc in its species.

Below you will find a short tutorial on how to use this pipeline in python. We will use the following example: We want to know which genes in the *Capitella* *teleta* genome are transcription factors (TFs).


## Required files
There are some files we will need in order to run this pipeline. The current version (Oct-2021) needs to gather some files manually, future versons will automatize this procedure.

### Query list
To run the pipeline, we need a list with human genes that represent all genes with the desired function in humans, in this case, a curated list with all human TFs. t do so, simply right a text/tsv/csv file with a single line per gene. Genes can be in any format we want. In this case we are using (the recommended) Ensembl GenID

### Eggnog database
We need to know for each of the genes in the "query" list, what are their respective orthogroups. This data is stored in the Eggnog database. We can use orthogroups at different levels. Let's say we want to relate Human and *Capitella* genes by both Metazoan and Bilaterian orthogroups. Then we need to goto the [Eggnog database](http://eggnog5.embl.de/#/app/home), go to Downloads, find out level of interest, in this case Metazoan and Bilaterian, one at a time, and download the (taxID)\_members.tsv.gz. De-compress it and save it in a folder.

### Lookup for conversions
The eggnog files we just downloaded show us all orthogroups witin a tax level and all **Proteins** that belongto that group. However, in most cases we won't have a list of **Proteins** but a list of **Ensembl Gen IDs** or something else. So we need to trannslate those Proteins in the eggnog files to whatever we have in our query list. To do so we need a lookup table that gives us for each protein in eggnog, its Ensembl Gen ID equivalent.

To do so we go to [Biomart](https://www.ensembl.org/biomart/martview/62f7cff64a6aaf9711bf7c0a3e52f7e7), and we select the species we want to translate, in this case, human. Then, in attributes we choose what are the conversions we want to make. We are hoosing **Protein stable** ID and **Gene stable ID**. Additionally, we want to know a more human readablename for each gene, so we also go to **EXTERNAL** and select HGNC symbol. Finally we click on Result in the top left of the screen and after it loads we hit GO! to download in tsv format.

> You could run the pipeline without the translation, for which, at the moment (Oct-2021)
> there is no in-built functionality, however, you can fake this by creating
> a lookup table where there are two identical columns, both with Protein IDs

### emapper output
Finally we need to know to which orthogroups each of our *Capitella* target genes belong to, in order to match this with the Eggnog database information. 
To do so, we need a proteome (or genome) in FASTA format, in our case *Capitella teleta*'s proteome. We run it using [emapper](http://eggnog-mapper.embl.de/) with standard values. We save the output file in a folder as is.


## The pipeline

### Install package
To do so follow the README.md file in our GitLab repository

### Import data
First we need to load all data we mentioned in the [Required data](##required-files) section. We wil show now how to import all these datasets.

1. Eggnog datasets
For this we have a function called read_eggnog, that allows you to import as many tax levels as you decide. Simply add as arguments a path for every (de-compressed) eggnog file you have.  
```
>>> eggnog = read_eggnog('/path/to/Eggnog_Bilateria(33213)_members.tsv', '/path/to/Eggnog_Metazoa(33208)_members.tsv')
``` 
2. lookup table
Simply import it using the pandas function read_csv. Don't forget to add whatever separator you used in *sep* argument ("\t" for tsv files, "," for csv, etc.)
```
>>> lookup = pandas.read_csv('/path/to/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')
```

3. List of Human transcription factors
Again, just a text/tsv/csv file with a single column containing all your curated genes. In our case, human transcription factors. Put a header in the file.
```
>>> query = pd.read_csv("/path/to/TFs_Ensembl_v_1.01.txt", sep = "\t")
```

4. emapper output
Target genes emaper output exactly as it comes from emapper. The argument *skiprows=4* removes the top 4 row of the file, containing metadata information. Do not remove this argument
```
>>> emapper = pd.read_csv("/path/to/Capitella_teleta_Prot_emapper_annotations.txt", skiprows=4, sep="\t")
```

### Run pipeline
First we need to create a table with all our human transcription factors and their respective orthogroups at all the taxonomic levels chosen. This can be breaken in two parts: Translating the eggnog database from protein IDs to gene IDs and then removing from there all proteins that are not in our query (list of human TFs).

First we translate the eggnog datasets with egg_translate() like this
```
translated_eggnog = egg_translate(eggnog, lookup)
```

Then we can run a quality control script to tell us wether any of the genes in the eggnog database were not translated
```
translated_QC(translated_eggnog)
```

If there were any non-translated genes, then the lookup table didnot catch all eggnog proteins, wether because they are outdated or anything else. 
We can check if any of the non-translated genes were not translated. If you find any of yur Curated genes is outputted by this function
you should go check in the eggnog database wether this isbecause your gene wasn't in the database at all (in which case, nothing should be done) or your gene was in the eggnog database but was not translated, in which case the lookup should be updated to account for that otherwise lost gene.
```
translated_QC(translated_eggnog, query)
```

Back to the pipeline, once we have checked the translation was optimal we can proceed by removing from the translated output any gene that is not a transcription factor. We use the following function. In merge_on we put the name of the column of translated_eggnog where we want to do the merge.
```
query_orthogroups = merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
```

This outputs a table with all your TFs, and their respective orthorgoups at the levels chosen like this
|Gene       |Orthogroup@33208|Orthogroup@33298|HGNC        |ProteinID      |
| ------ | ------ | ------ | ------ | ------ |
|Query Gene1|FGFV@33208      |BVGC@33298      |Gene symbol1|Target Protein1|
|Query Gene2|TGOK@33208      |CACF@33298      |Gene symbol2|Target Protein2|

> The keep_conversions argument allows the user to remove any other gene conversions 
> introduced by the lookup table. In our case it would remove the HGNC and ProteinID columns

The Orthogroup@ columns show you the orthogroups at the different taxonomical levels. the level is denoted by the @(taxID)

If we combine that with the emapper output we can know which of the *Capitella* genes (emapper) share an orthogroup at Metazoan or Bilaterian level with a Human TF. We can do this with the following code:
```
annotated_genes = emapper_annotation(emapper, query_orthogroups, keep_all_targets= False)
```

annotated_genes contains a table with all *Capitella* genes that are putative TFs, based on your list of Human Tfs and the taxonomic evels chosen at Eggnog. You can save this file as a tab separated file with the following code:
```
annotated_genes.to_csv("/g/arendt/Javier/Python/geneannotator/tests/Ortho_method_Capitella_TFs.tsv", sep = "\t")
```
