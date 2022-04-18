#####################################
#### Use emapper output to find out which 
#### of your emapper genes have a certain GO:Term
#### and therefore are a particular gene module
#####################################


import pandas as pd
import numpy as np
pd.options.display.max_rows = 999
pd.options.display.max_columns = 999



def GOTerms_annotation(emapper, GOterm, extra_columns=False, keep_all_columns = False):
	'''
	This function takes in emapper results and annotates genes as TFs if they have a specified GO:Term (for TFs is GO:0003700)
	
	Attributes
	----------
	emapper: string
		Path to emapper output file. emapper file should not be modified
	GOterm: string
		GO:term you want to search on the mapper output
	'''

	emapper = pd.read_csv(emapper, skiprows=4, sep="\t")

	# What columns should we keep in final output?
	columns = ["#query"]
	if isinstance(extra_columns, list) and keep_all_columns == False:
		columns.extend(extra_columns)

	elif keep_all_columns:
		columns = emapper.columns.values

	elif extra_columns is False:
		pass
		#print("No emapper extra_columns added to final output")

	else: exit("extra_column is not a list, please input a list")



	# Contains extra_columns any wrong name?
	if len(emapper.columns.intersection(columns)) != len(columns):
		exit("You inserted a wrong column name in extra_columns")

	annotation = emapper["GOs"].str.contains(GOterm)
	annotation.fillna(False, inplace = True) ## Last 3 rows of emapper are statistics
	annotation = emapper[emapper.columns.intersection(columns)][annotation]
	
	# Re-format dataframe
	annotation = pd.DataFrame(data = annotation)
	annotation.reset_index(inplace = True)
	annotation.drop(columns = ["index"], inplace = True)
	
	return pd.DataFrame(data = annotation)




##### Scriptexample usage
"""
pd.set_option('display.max_columns', None)
print(GOTerms_annotation("/g/arendt/Javier/Python/TF_annot_methods/Capitella_teleta/Data/Capitella_teleta_Prot_emapper_annotations.txt", ["EC", "GOs"], GOterm="GO:0003700"))
"""