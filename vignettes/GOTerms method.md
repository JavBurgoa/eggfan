# **GO:Terms method**
This method annotates an emapper output based on the GO:Terms specified in this very output. It simply tells you which genes in your emapper output have a particular GO:Term of your choice. This allows you to know which genes belong to a particular module/family (like for example which genes are transcription factors) as long as you know the GO:Term that identifies that module. More information about GO:Terms can be found in the [Gene ontology resource website](http://geneontology.org/)

For demonstration purposes we will keep our example of:  Which genes in the *Capitella* *teleta* genome are transcription factors (TFs)?.

&nbsp;
&nbsp;

## **Required files**
### Emapper output
We need to know to whihc GO:Terms our target genes belong to,  
To do so, we need a proteome (or genome) in FASTA format, in our case *Capitella teleta*'s proteome. We run it using [emapper](http://eggnog-mapper.embl.de/) with standard values. We save the output file in a folder as is. If you just want to try the algorithm you can use the file in /test/data/Capitella_emapper_redux.txt .

### GO:Term
This is not a file, it is just the GO identifier that encapsulates the gene module we want to annotate. In this case the GO:Term for transcription factors is `GO:0003700`:



&nbsp;
&nbsp;


## **The pipeline**

### Python
### Install package
To do so follow the README.md file in our GitLab repository

### Import modules
For this pipeline you don't need to import any file. Just open your script, and import the following packages
```
>>> from eggfan import goterms
```

### Run pipeline
Now, to run the pipeline, as long as you know your GO:Terms and the path to your emapper output you should be able to run
```
>>> cap_annotated = goterms.GOTerms_annotation("path/to/capitella_emapper.txt", GOterm="GO:0003700")
>>> print(cap_annotated)
```

This should output a list with all the genes/proteins in the ampapper file that have that GO:Term.
there are some optional parameters one can add. `extra_columns` flag allows you to input a list with all the extra columns in the emapper output you want to include other than the genes column.

```
>>> cap_annotated = GOTerms_annotation("path/to/capitella_emapper.txt", GOterm="GO:0003700", extra_columns=["GOs", "ECs"])
>>> print(cap_annotated)
```
Finally if you wish to keep all columns from emapper output you can use the flag `keep_all_columns = True`. 

&nbsp;
&nbsp;

### **Command line**
You can see all options in the command line orthogroup function simply with
```
python src/eggfan/goterms_argparse.py -h
```


1. **-em**
Path to emapper output run with your target proteins you want to annotate [Emapper output](###-emapper-output).
If you just want to try the algorithm you can use the file in /test/data/Capitella_emapper_redux.txt .
2. **-g**
Goterm code you want to search in the emapper output.

#### Run pipeline
Simply run the pipeline like so:

```
python eggfan/goterms_argparse.py -e /test/data/Capitella_emapper_redux.txt -g "GO:0003700"
```
This outputs a list with all the genes/proteins in yout emapper that contain the specified GO:Term.
There are some extra flags that modify the output. `-x` allows you to input emapper column names you would like to see in the output. For example if you also want to see the other GO:terms and eggnog Orthogroups yo can run:

```
python eggfan/goterms_argparse.py -e "/test/data/Capitella_emapper_redux.txt" -g "GO:0003700" -x"GOs" -x"eggNOG_OGs"
```
Additionally with the `--keep_all_columns` flag you can keep all emapper output columns.