###############
#### This script uses the orthology table from the phylome
#### A resource containing all genes of a particualr species and speicfying all genes in any other species that are orthologous to it
#### Based on this and knowing which human genes are belong to a certain module (TFs, contractile genes etc), we can say that if a gene has as ortholog a human TF, 
#### then the gene it also belongs to that "module"
###############

import pandas as pd
import numpy as np
import os
import utils
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 500)

def initial_lookup(path):
	"""
	Based on all the genes from the inputted orthology table(s), make a lookup that translates all of them from UNIrpot ID to ENSEMBL ID and HGNC

	Attributes
	----------
	path: string.
		Path to phylome orthology table(s)
	"""
	genes = utils.human_genes_string(path)

	uni_ens_lookup = utils.uniprot_request(genes, 'ID', "ENSEMBL_ID")
	ens_HGNC = utils.uniprot_request(genes, 'ID', 'GENECARDS_ID')

	uni_ens_lookup = utils.format_uniprot_output(uni_ens_lookup, ["UniProtKB", "ENSEMBL_ID"])
	ens_HGNC = utils.format_uniprot_output(ens_HGNC, ["UniProtKB", "HGNC"])

	lookup = uni_ens_lookup.merge(ens_HGNC, how="outer", left_on="UniProtKB", right_on="UniProtKB")

	return lookup, genes




def check_lost_genes(genes, lookup):
	"""
	Make a table with all Human uniprotIDs that were translated either to ENSEMBL_ID or HGNC or nothing, but not to the both of them (so, the ones where some translation is missing)
	
	Attributes
	----------
	genes: string
		All human Unirpot IDs to be translated, this is the IDs in th orthology tables
	lookup: pandas dataframe
		lookup table with UnirpotKBs, ENSEMBL_ID and HGNC as columns with some missing translations to be completed
	"""

	human_genes = genes.split(" ")
	
	translated = list(lookup["UniProtKB"])
	lost_genes = pd.DataFrame(human_genes)[~np.isin(human_genes, translated)]
	lost_genes.columns= ["UniProtKB"] # genes that wewre not tranlated to anything

	no_genID = lookup[lookup["ENSEMBL_ID"].isna()]
	no_HGNC = lookup[lookup["HGNC"].isna()]

	lost_genes = lost_genes.merge(no_genID, how="outer", left_on="UniProtKB", right_on="UniProtKB")
	lost_genes = lost_genes.merge(no_HGNC, how="outer", left_on="UniProtKB", right_on="UniProtKB")
	
	lost_genes.fillna("", inplace=True)
	lost_genes["Translation"] = lost_genes.iloc[:, 1] + lost_genes.iloc[:, 3] + lost_genes.iloc[:, 2]
	lost_genes.drop(columns=["ENSEMBL_ID_x", "HGNC_x", "ENSEMBL_ID_y", "HGNC_y"], inplace = True)

	return lost_genes


def correct_uniprot_translation(up):
	"""
	Sometimes, you ask Uniprot to give you HGNC and it gives you ENSEMBL, so, find where that happens in the original lookup and put those errors inplace (also, when Uniprot does this, it doesnt find the ENsemble translation)
	
	Attributes
	----------
	lookup: pandas dataframe
		The same lookup table as before but the updated version from translate_from_HGNC(). Could be used with the non updated version aswell.
	"""


	return lookup


def translate_from_HGNC(lost_genes, lookup):
	"""
	For those human uniprots that only got HGNC translation, take that HGNC,  translate it to ENSEMBL_ID and update the lookup with that information
	
	Attributes
	----------
	lost_genes: pandas dataframe
		product of check_lost_genes(). A dataframe with all Human UnirptoIDs from the orthology tables taht were not fully translated to ENSID or HGNC and the respective translation they were given
	lookup: pandas dataframe
		lookup table with UniprotKBs, ENSEMBL_ID and HGNC as columns with some missing translations to be completed
	"""

	# For those human uniprots that only got HGNC translation, take that HGNC and translate it to ENSEMBL_IDs
	lost_genes = lost_genes.replace(r'^\s*$', np.nan, regex=True) # replace empty strings with NAs, for later dropna()
	lost_genes = lost_genes.dropna()
	table = {}
	for i in lost_genes.index.values:

		HGNC = lost_genes.loc[i, "Translation"]
		Uniprot = lost_genes.loc[i, "UniProtKB"]

		if not HGNC.startswith("ENSG0000"):
			
			ENSG = utils.HGNC_request(gene = HGNC)
			table[i] = [Uniprot, ENSG, HGNC]
			print(table[i])

	table = pd.DataFrame.from_dict(table, orient = 'index').dropna()
	table.columns = lookup.columns
	Updated = lookup.merge(table, how = "left", left_on = "UniProtKB", right_on="UniProtKB")	

	# Sometimes, you ask Uniprot to give you HGNC and it gives you ENSEMBL, so, find where that happens in the original lookup and put those errors inplace (also, when Uniprot does this, it doesnt find the ENsemble translation)
	rows_without_genid = Updated[Updated["ENSEMBL_ID_x"].isna()].index.values
	for i in rows_without_genid:
		# sometimes, you ask Uniprot to give you HGNC and it gives you ENSEMBL, so put them inplace (also, when Uniprot does this, it doesnt find the ENsemble trnaltion)
		if Updated.loc[i, "HGNC_x"].startswith("ENSG0000"):
				Updated.loc[i, "ENSEMBL_ID_x"] = Updated.loc[i, "HGNC_x"]
				Updated.loc[i, "HGNC_x"] = ""
		# Add geneIDs found in HGNC
		Updated.loc[i, "ENSEMBL_ID_x"] = Updated.loc[i, "ENSEMBL_ID_y"]

	Updated = Updated.drop(columns = ["ENSEMBL_ID_y", "HGNC_y"])
	Updated.columns = lookup.columns
	Updated.drop_duplicates(inplace = True)

	return Updated








# Tranlsated method

###########################
#### Make Lookup table ####
def make_lookup(path, out = False):
	"""
	Gets all uniprot IDs from the specified orthology tables and makes a lookup table that translates them to whatever you desire. Default Ensembl IDs.

	Attributes
	----------
	path: string.
		Path to phylome orthology table. Path can be a file path or a path to a folder. In the latter case it will run the pipeline for the whole folder
	out: string
		path to file where you want to save the lookup. Default, not saving anything
	"""
	# make lookup table and a list with all Unirpot IDs to be translated
	lookup, uniprots = initial_lookup(path)

	# make list with all genes not fully translated (either ENSEMBL_ID or HGNC were not retrieved)
	lost_genes = check_lost_genes(uniprots, lookup)

	## Update table by translating HGNCs to ENSEMBLIDs when possible
	updated = translate_from_HGNC(lost_genes, lookup) # this is very slow beacuase genecards only allows one gene at a time to translate

	if out is not False:
		updated.to_csv(out, sep="\t", index=False)

	return updated




# Make translated orthology tables
def translate_orthologies(path, lookup, out = False):
	"""
	Takes in one or several phylome orthology tables and translates their human UniprotIDs to ENSEMBL and HGNC, adding an extra column on each of the orthology tables inputed. Output is a list with a dataframe per orthology table
	path: string.
		Path to phylome orthology table. Path can be a file path or a path to a folder. In the latter case it will run the pipeline for the whole folder
	lookup: pandas dataframe
		output from make_lookup(). A lookup table with three columns. ENSEMBL_ID, HGNC and UniProtKB, with the translations of human genes in each of those ID types.
	out: string (optional)
		path to DIRECTORY where you want the file(s) to be saved in case you are using various files, in shih¡ch case they should have the default name taxID_orthogroup.tsv . They will be given a slightly different name than the original by default, adding the suffix "_human_". If you just have one file you can specify the output name in the path
	"""
	orthology_tables = utils.directory_or_file(path)
	
	tables = []
	for fullpath in orthology_tables:
		# Import
		orthoTable = pd.read_csv(fullpath, index_col=False, skiprows=[i for i in range(1,13)], sep = "\t")
		
		# Format orthotables and lookup
		orthoTable = orthoTable[orthoTable["target_species"] == "Homo sapiens"]
		orthoTable["ENSEMBL_ID"] = ""

		# Translate
		orthoTable = orthoTable.fillna("")
		translated_orthoTable = utils.translate_uniprots(orthoTable, lookup)

		# Format
		translated_orthoTable["ENSEMBL_ID"] = translated_orthoTable["ENSEMBL_ID"].replace("^,", "", regex=True) # the pipeline added an extra comma in the beginning by default
		translated_orthoTable["ENSEMBL_ID"] = translated_orthoTable["ENSEMBL_ID"].replace(",,", ",-,", regex = True).replace(",,", ",-,", regex = True).replace("^,", "-,", regex = True).replace(",$", ",-", regex = True) # add dashes where missing genes (a relpace is repeated on purpose)
		
		tables.append(translated_orthoTable)
	
	# Save
	if os.path.isdir(out):
		for i in range(len(orthology_tables)):
			filename = os.path.basename(orthology_tables[i])
			file = out + filename.replace("_orthologs.tsv", "_human_orthologs.tsv")
			tables[i].to_csv(file, index = False, sep = "\t")

	elif isinstance(out, str):
		tables[0].to_csv(out, index = False, sep = "\t")
	

	return tables





def read_translated_tables(translated_orthologies):
	"""
	Make a list of pandas dataframes if the input is a directory
	Keep as is if input is a list of dataframes
	"""
	if isinstance(translated_orthologies, str):
		orthology_tables = []
		for file in os.listdir(translated_orthologies):
			table = pd.read_csv(translated_orthologies + file, sep = "\t")
			orthology_tables.append(table)

	elif isinstance(translated_orthologies, list):
		orthology_tables = translated_orthologies
		del translated_orthologies
		
	else: exit("translated_orthologies input is not a list nor a directory path")

	return orthology_tables





def subset_query_orthologs_and_position(orthoTable, human_query):
	"""
	Ths function takes the orthology tables (translated) and subsets them to include only the genes that have as orthologs ghuman genes in our query.
	Additionally, it creates a separate table specifying in which position in the orthology those genes were. Witth this I mean:

	If you have:
	|    Seed   |           ortholog        |
	1 | 456.GeneA | HumanA, HumanA7           |
	2 | 768.GeneB | HumanB4, HumanB3, HumanB6 |
	
	and Human A7 and HumanB4 are in our query, this function would make a table of their position like this:
	
	  | "GenID" | "position_from_0" | "number_of_IDs" |
	1 | HumanA7 | 1                 | 2				  |
	2 | HumanB4 | 0                 | 3               |

	Then this table will be used to add extra columns saying whihc of all the orthologs in the original phylome tables are actuaolly in our query. And they will match the position in the original table like so:
	
	  |    Seed   |           ortholog        | frrom_query |
	1 | 456.GeneA | HumanA, HumanA7           | -,HumanA7   |
	2 | 768.GeneB | HumanB4, HumanB3, HumanB6 | HumanB4,-,- |


	Attributes
	----------
	orthoTable: pandas dataframe
		Translated orthology table product of translate_orthologies()

	human_query: pandas dataframe
		Containing a single column with all ENSEMBL gene IDs. They represent a module/family in humans that you want to identify in your target species. 
	"""
	
	finalorthotable = pd.DataFrame(columns = {"##Seed_(co-)orthologs":[], "type":[], "ENSEMBL_ID":[], "orthologs": [], "GeneName_target":[], "ENSEMBL_query-only":[], "GeneName_target_query-only":[]})
	query_position = pd.DataFrame({"GenID":[], "HGNC":[], "position_from_0":[], "number_of_IDs":[]})

	# Create two dataframes: One with all phylome rows that contain TFs and one with position information on each one of the TFs found in the orthotable
	for gene in human_query.iloc[:, 0].unique():

		# Get all orthology tables rows that contain |a gene in the query| as ortholog
		condition = orthoTable["ENSEMBL_ID"].str.contains(gene)
		condition.fillna(value=False, inplace = True) # !!!Can be elimnated now that we have "-"??. those are one to one orthologs that were not translated
		rows_perGene = orthoTable[["##Seed_(co-)orthologs", "type", "orthologs", "GeneName_target", "ENSEMBL_ID"]][condition]

		# Make tables
		if len(rows_perGene) > 0:
			# add final subset per human ortholog
			finalorthotable = finalorthotable.append(rows_perGene)

			# Create table with which gene was found in which position in which row (in a row with various orthologs separated by commas, which one is the one that matches the query)
			position = utils.find_position(rows_perGene, gene, column = "ENSEMBL_ID")
			table = utils.query_position_table(rows_perGene, position, gene)
			query_position = query_position.append(table)

	# Final formatting
	finalorthotable = finalorthotable[~finalorthotable.index.duplicated(keep='first')]
	query_position = query_position.rename_axis('idx').sort_values(by = ['idx', 'position_from_0'])

	return finalorthotable, query_position



def HGNC_subset_query_orthologs_and_position(orthoTable, human_query):
	"""
	This function is sister to subset_query_orthologs_and_position(). Same as that one takes the orthology tables (translated) and subsets them to include only the genes that have as orthologs ghuman genes in our query.
	Additionally, it creates a separate table specifying in which position in the orthology those genes were. More infor in subset_query_orthologs_and_position() docs
	This function has certain specificities to it that I think ould make the merging of these two functions a bit messy, despite how simmilar they are.

	Attributes
	----------
	orthoTable: pandas dataframe
		Orthology tables without any adulterations, straight from phylome output

	human_query: pandas dataframe
		Containing a single column with all HGNC (Genecards) gene IDs. They represent a module/family in humans that you want to identify in your target species. 
	"""
	symbol_col_name = "GeneName_target"
	finalorthotable = pd.DataFrame(columns = {"##Seed_(co-)orthologs":[], "type":[], "GeneName_target":[], "GeneName_target_query-only":[]})
	query_position = pd.DataFrame({"HGNC":[], "position_from_0":[], "number_of_IDs":[]})
	
	for gene in human_query.iloc[:, 0].unique():
		condition = orthoTable[symbol_col_name].str.contains("(?:^" + gene +"$|^" + gene + ",|," + gene +",|," + gene + "$)", regex = True)
		condition.fillna(value=False, inplace = True) # fill in one to one orthologs that were not trnalated
		rows_perGene = orthoTable[["##Seed_(co-)orthologs", "type",  "GeneName_target"]][condition]

		if len(rows_perGene) > 0:
			# add final subset per human ortholog
			finalorthotable = finalorthotable.append(rows_perGene)

			# Create table with which gene was found in which position in which row
			position = utils.find_position(rows_perGene, gene, column = symbol_col_name, HGNC = True)
			table = utils.query_position_table(rows_perGene, position, gene, HGNC = True)
			query_position = query_position.append(table)
	
	finalorthotable = finalorthotable[~finalorthotable.index.duplicated(keep='first')]
	query_position = query_position.rename_axis('idx').sort_values(by = ['idx', 'position_from_0'])

	return finalorthotable, query_position



def add_queryonly_columns(finalorthotable, query_position, HGNC = False):
	"""
	Add to the final tables columns specifying which of the genes in the orothologs column appears in the human query

	Attributes
	----------
	finalorthotable: dataframe
		translatedorthology tables from phylome wcontaining only genes with humna orthologs in our query and two empy columns to be filled by this function
	query_position
		table specific for each orthology table specifying which of the genes in the orothologs column appears in the human query and their position in the set. This is explained better in the docs from subset_query_orthologs_and_position()
	HGNC: Boolean
		Whether you want the version for the norml pipelien of the HGNC version
	"""

	index = finalorthotable.index.values
	for i in index:
		rows = query_position[query_position.index == i] # make a dataframe with all the TFs found in one row (i) in finalorthotable

		rows.set_index("position_from_0", inplace=True)
		new_index = list(range(int(rows.number_of_IDs.iloc[0])))

		rows = rows.reindex(new_index, fill_value="-") # If a gene was not a Tf substitue it with "-"
		rows.index= list(range(len(rows.index.values)))

		if not HGNC:
			finalorthotable["ENSEMBL_query-only"][finalorthotable.index == i] = ','.join(list(rows["GenID"]))
		finalorthotable["GeneName_target_query-only"][finalorthotable.index == i] = ','.join(list(rows["HGNC"]))

	# remove this after adding drop_duplicates in the very bginning.:
	finalorthotable = finalorthotable.drop_duplicates() # It would be better to put this drop duplicates in the beginning, for each translated table, to make things a bit faster.

	return finalorthotable




def find_query_orthologs(query_path, translated_orthologies):
	"""
	Find in phylome (all genes/proteins of a target species) which genes/proteins have as orthologs any gene/protein in your human query. Make a table out of it.

	Attributes
	----------
	query_path: string
		Path to list of human ENSEMBL IDs that represent the gene module/family you want to search in phylome species.
	translated_orthologies: list or path
		Phylome orthology tables with human orthologs Unitrots translated to ENsemblIDs. This is, the product of translate_orthologies(). You can input a path to a folder containing all of those tranlslate orthologies or a lists object full of pandas dataframes.

	"""

	## Import data
	human_query = pd.read_csv(query_path)
	orthology_tables = read_translated_tables(translated_orthologies)

	# Subset and add columns
	tables = []
	for orthoTable in orthology_tables:

		finalorthotable, query_position = subset_query_orthologs_and_position(orthoTable, human_query)
		
		finalorthotable = add_queryonly_columns(finalorthotable, query_position)

		tables.append(finalorthotable)


	return tables

def save_annotated(annotated_tables, directory, suffix = "_annotated_orthology"):
	"""
	Saves a list of dataframes into separate dataframes with specific names

	Attributes
	----------
	annotated_tables: list
		list containing all dataframes to be saved
	directory: string
		path to directory where you want to save the files
	suffix: string
		name of output file will be <taxID><suffix>.tsv . Default "_annotated_orthology"
	"""
	if directory[-1] != "/":
		directory = directory + "/"

	for table in annotated_tables:
		taxID = table.iat[0, 0].split(".")[0]
		file = directory + taxID + suffix + ".tsv"

		table.to_csv(file, index = False, sep = "\t")



############################
###### HGNC method #########

# Make final tables with GeneID(s) species | orthology type | all human Ensembl orthologs | All HGNCs | TF EnsemblIDs | TF HGNCs
def annotate_orthology_HGNC_method(query_path, orthology_tables_path):
	"""
	orthology_tables_path: string
		path to folder containing orthology tables you want to annotate. Alternatively you can input a path to a single file
	"""
	## Import data
	human_query = pd.read_csv(query_path)
	orthology_tables_path = utils.directory_or_file(orthology_tables_path)

	tables = []
	for file in orthology_tables_path:
		# Import
		orthoTable = pd.read_csv(file, index_col=False, skiprows=[i for i in range(1,13)], sep = "\t")

		# Format
		orthoTable = orthoTable[orthoTable["target_species"] == "Homo sapiens"] ## Added for the HGNC pipeline

		# Subset dataframes and locate query genes
		finalorthotable, query_position = HGNC_subset_query_orthologs_and_position(orthoTable, human_query)

		#### HGNC exclusive part ####
		# Remove duplicates, because thep ipeline somehow duplicates the genes found in position 0. But account for the index, in case sma ortholog found in same position but itn different line
		query_position['index'] = query_position.index
		query_position.drop_duplicates(inplace=True)
		del query_position['index']
		##############################

		finalorthotable = add_queryonly_columns(finalorthotable, query_position, HGNC = True)

		tables.append(finalorthotable)

	return tables



### Script
"""
#lookup = make_lookup("/g/arendt/data/phylomeV2/orthology_tables_nocollapse/5759_orthologs.tsv")
lookup = pd.read_csv("/g/arendt/Javier/Python/geneannotator/tests/lookup_5759.txt", sep = "\t", keep_default_na=False)
translated_orthologies = translate_orthologies("/g/arendt/Javier/Python/geneannotator/tests/translated_orthotables/", lookup)
annotated_tables = find_query_orthologs("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", translated_orthologies)
#print(annotated_tables)
# Save
save_annotated(annotated_tables, "/g/arendt/Javier/Python/geneannotator/tests/")
"""


### HGNC version
"""
annotated_tables = annotate_orthology_HGNC_method("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TF_names_v_1.01.txt", "/g/arendt/Javier/Python/geneannotator/tests/translated_orthotables/")
save_annotated(annotated_tables, "/g/arendt/Javier/Python/geneannotator/tests/")
"""