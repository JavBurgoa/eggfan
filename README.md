# Eggfan
![Venn diagram](eggfan.png "The three pipelines complementarity")
<br></br>


EggFan (**Egg**nog based **F**unctional **an**notator) consists of three separate but complementary pipelines that annotate the function of any gene/protein of your interest, using some publicly resources. These resources are [Egggnog](http://eggnog5.embl.de/#/app/home), [emapper](http://eggnog-mapper.embl.de/) and the Phylome, all resources made and mantained by the Jaime Huerta Cepas lab in the CBGP unit of Universidad Politécnica de Madrid.

If you want to know more about each pipeline, check out the vignettes in the github repository.

# Installation

## pip install (**Not ready!**)
The easiest way is just running on your terminal
```
pip install eggfan
```

## Manual installation

If you don't already have one, create a Conda / Miniconda environment:
```
conda create -n eggfan -c conda-forge python=3.7 pip pandas numpy
```

Install eggfan
```
git clone https://git.embl.de/burgoa/geneannotator.git 
cd eggfan
pip install .
```
