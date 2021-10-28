import pandas as pd
import numpy as np




def read_eggnog(*paths):
    """
    Read and format one or many eggnog datasets for posterior uses by further functions. Formatting consists of adding column names
    ...

    Attributes
    ----------
    *paths: string(s)
        Absolute or relative path(s) to eggnog datasets as downloaded from eggnog (Go to http://eggnog5.embl.de/#/app/downloads, click on the taxonimic level you are interested in, and download the file with suffix "members.tsv.gz")
    """
    eggnogs = []
    for path in paths:
        eggnog = pd.read_csv(path,
                         sep='\t', 
                         names=["X", "Orthogroup", "N_Prots", "N_Spec", "Protein stable ID", "SpeciesID"]
                         )
        eggnogs.append(eggnog)


    if len(eggnogs) > 1:
        return eggnogs
    else:
        return eggnogs[0]


# This fucntion is to make a future lookup table using APIs
def eggnog_prots_extract(eggnog, taxID):
    """
    Find all proteins from a species in an eggnog dataset
    ...

    Attributes
    ----------
    eggnog : pandas dataframe
        eggnog dataset as downloaded from eggnog and formatted by read_eggnog()
    taxID: string
        NCBI tax ID of the species you want to retrieve proteinIDs from. default = human
    """

    all_proteins = ""

    eggnog = eggnog.ProtID.str.split(',')
    for array in eggnog:
        prot = [x for x in array if x.startswith(taxID)]
        if len(prot) > 0:
            prot = ",".join(prot)
            all_proteins = all_proteins + "," + prot

    # format IDs
    all_proteins = all_proteins.replace(taxID + ".", "")

    # remove duplicates
    all_proteins = all_proteins.split(",")
    all_proteins = set(all_proteins)
    all_proteins = ','.join(all_proteins)

    return all_proteins

def query_table(dataset, match_column, taxID, data_origin):
	"""
    Takes in a table with two columns, one of them with unique identifiers and the other with multiple identifiers separated y commas as:
    | protein | Orthogroup                        |
    | ENSP01  | GHFC@metaz, GBHF@metaz, GFCB@Opis |
    | ENSP02  | FVHG@Bilate, HJHG@Opis, HGNJ@Phot |

    And an id that identifies some of the multiple identifiers (Orthogroup columns in this example). It reduces the table size by picking only row that contain the id
    and removes all identifiers that do not contain the id. thus for the id = metaz we would get:

    | protein | Orthogroup            |
    | ENSP01  | GHFC@metaz, GBHF@metaz|

    Attributes
    ----------
    dataset: pandas dataframe
    	table with two columns as explained in the documentation above. It is relevant that one of the columns contains unique identifiers and the other multiple.
    match_column: string
    	name of the column with the multiple identifiers where we want to match the taxID
    taxID: string
    	substring, present in the match_column column. All identifiers without the TaxID will be eliminated
    data_origin: string
    	can take values "emapper" or "eggnog", depending on whether we want to do this on a dataset downloaded from eggnog or the output of emapper. emapper input is usually two columns, one
    	with protein IDs from the target species ("#query") column and one with all orthorgoups from that protein separated by ",". eggnog is simmilar, one column with one orthorgoup per row and another one
    	with several proteins, separated one from another by ",".
    """
	draged_column = [colname for colname in dataset.columns.values if colname != match_column][0]
	if isinstance(dataset.loc[0, match_column], str): # This is a patch. when using this function in a loop like emaper_annotation() the original emapper object gets modified by this line the first time. If you do it again you get NAs. I don't know why this happens
		dataset.loc[:, match_column] = dataset.loc[:, match_column].str.split(',')
	out = pd.DataFrame()

	for i in dataset.index.values:
		if data_origin == "emapper": # not very efficient evaluating for every row
			matched_element = [element for element in dataset.loc[i, match_column] if taxID in element]
		elif data_origin == "eggnog":
			matched_element = [element for element in dataset.loc[i, match_column] if element.startswith(taxID)]

		if len(matched_element) > 0:
			matched_element = ",".join(matched_element)
			matched_element = pd.DataFrame(data = {draged_column:[dataset.loc[i, draged_column]], match_column : [matched_element]})
			out = out.append(matched_element, ignore_index = True)

	return(out)



def eggnog_orthoprot_table(eggnog, taxID, explode = True, remove_taxid = True):
    """
    Find all proteins from a species and respective orthogroup in an eggnog dataset
    ...

    Attributes
    ----------
    eggnog : pandas dataframe
        eggnog dataset as downloaded from eggnog and formatted by read_eggnog()
    taxID: string
        NCBI tax ID of the species you want to retrieve proteinIDs from. default = human
    explode: boolean
        If true it will make a single row per protein ID. Otherwise all proteins that were in the same row, will remain in the same row separated by ","
    remove_taxid: boolean
        If true it will remove taxID from the prot IDs. Otherwise it will leave them with their protID
    """

    prot_column = "Protein stable ID"

    prot_ortho = query_table(eggnog, match_column = prot_column, taxID = "9606", data_origin = "eggnog")
    prot_ortho.drop_duplicates(inplace = True)

    # Return
    if remove_taxid:
        prot_ortho.loc[:, prot_column] = prot_ortho.loc[:, prot_column].str.replace(taxID + ".", "", regex=False)

    if explode:
        prot_ortho.loc[:, prot_column] = prot_ortho.loc[:, prot_column].str.split(",")
        prot_ortho = prot_ortho.explode(prot_column)

    return prot_ortho


def egg_translate(eggnog, lookup, taxID = "9606"):
    """
    Takes in one or several Eggnog database raw datasets , finds each of the proteins in each eggnog dataset and creates a table specifying for each protein their respective Ensembl gene ID, HGNC symbol (HGNC and Ensembl gen ID may not be "translated", Na for empty values) and orthogroup.
    ...

    Attributes
    ----------
    eggnog : list or pandas dataframe
        list containing one or more eggnog datasets in pandas dataframe format or a single pandas dataframe. Supposed to be direct output from read_eggnog.
    lookup : pandas.dataframe
        DataFrame containing three columns: "HGNC symbol", "Gene stable ID", "Protein stable ID". Each column contain strings with ID conversions from HGNC to Ensembl GenID to Ensembl protein ID
    
    Output
    ------
    - If input is a single eggnog dataset the output is a single table 
    - If input is a list with multiple eggnog datasets then the output is a single table with an extra column per dataset specifying orthogroups from the differnt datasets inputed
    """

    prot_column = "Protein stable ID"

    # Data preparation
    if not isinstance(eggnog, list):
        eggnog = [eggnog]

    tax_levels = []
    for df in list(range(len(eggnog))):
        tax = "@" + str(eggnog[df].iat[0,0])
        tax_levels.append(tax) # For later use in naming orthogroup columns

        eggnog[df] = eggnog[df][eggnog[df].SpeciesID.str.contains(taxID)] # only rows with human prots stay.
        eggnog[df] = eggnog[df].loc[:, [prot_column, "Orthogroup"]]

    lookup.dropna(subset=[prot_column, "Gene stable ID"], inplace=True)

    ## Make translated table(s)
    dfs = []
    for egg in eggnog:
        egg_prots = eggnog_orthoprot_table(eggnog = egg, taxID = taxID)
        egg_prots = egg_prots.merge(lookup, how = "left", on = prot_column)
        dfs.append(egg_prots)

    # Save
    out = dfs[0]
    if len(dfs) > 1:
        for df in dfs[1:]:
            df = df.drop(columns = ["HGNC symbol", "Gene stable ID"])
            out = out.merge(df, on = prot_column, suffixes = tax_levels)
    return out



def translated_QC(translated_eggnog, query = False, match_column = "Gene stable ID"):
    ortho_cols = [colname for colname in translated_eggnog.columns.values if colname.startswith("Orthogroup")]

    translated_eggnog = translated_eggnog.drop(columns = ortho_cols)
    lost = translated_eggnog.loc[translated_eggnog[match_column].isna(), :]
    lost_prots = lost["Protein stable ID"]
    print("Proteins that were in eggnog but were not translated:")
    print(lost_prots)

    if query is not False:
        query.columns = [match_column]
        query = np.array(query[match_column])
        print("Genes from your query that are not translated or were not in eggnog database:")
        lost_genes = np.array(lost[match_column])
        print(query[np.isin(query, lost_genes)])





def merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions = False):
    """
    Subset your genes of interest (query) from the translated_eggnog table. Outputs all orthogroups for your genes of interest + (if keep_conversions = True) conversions to other symbols
    ...

    Attributes
    ----------
    translated_eggnog: pandas dataframe
        Product of egg_translate with one or more eggnog datasets.
    query: pandas dataframe
        Dataframe containing a single row with Genes or proteins in ID format matching any column in translated_eggnog (Ensembl GenIDs, Ensebl Protein IDs, or HGNC symbols)
    merge_on: string
        String with the name of the column in lookup and traslated_eggnog you want to merge. Default: By Ensembl Gen IDs
    keep_conversions: boolean
        Wether to keep the rest of the translations of the genes (HGNC, Ensembl Gene ID and Ensembl Protein ID)
    """
    query.columns = [merge_on]
    query_with_orth = translated_eggnog.merge(query, how = "right", on = merge_on)
    ortho_cols = [colname for colname in query_with_orth.columns.values if colname.startswith("Orthogroup")]
    query_with_orth.dropna(subset=ortho_cols, inplace = True)

    if not keep_conversions:
        columns_keep = ortho_cols + [merge_on]
        query_with_orth = query_with_orth.loc[:, columns_keep]

    return query_with_orth

def format_quer_orth(query_orthogroups, ortho_cols):
	"""
	query_orthogroups has each orthogroup in a separate column. This script puts eveything in a single "Orthogroup" column, adding the @tax_ID to each orthogroup.

	Attributes
    ----------
    query_orthogroups: pandas dataframe
		Pandas dataframe with genes from query list, their conversions and orthogroups they belong to (in separate columns)
	ortho_cols: array
		Names of columns that contain the orthogroups 
	"""
	non_ortho_cols = [element for element in query_orthogroups.columns.values if element not in ortho_cols ]
	out = pd.DataFrame()

	for col in ortho_cols:
		# Make dataframe with only one Ortho column but all gene names and conversions
		non_ortho_cols = [element for element in query_orthogroups.columns.values if element not in ortho_cols ]
		subset = non_ortho_cols
		subset.append(col)
		new_query_orth = query_orthogroups.loc[:, subset]

		# Add @ taxID to orthogroups
		orthoname = col.replace("Orthogroup", "")
		new_query_orth.loc[:, col] = new_query_orth.loc[:, col] + orthoname

		new_query_orth.columns = ["Orthogroup" if element==col else element for element in subset] # we rename Orthogroup column to just "Orthorgoup" for future append
		out = out.append(new_query_orth)
	
	out = out.reset_index()
	
	return out

def format_query_targets(query_targets):
	group_by_col = "#query"
	non_query_cols = [colname for colname in query_targets.columns.values if not colname == group_by_col]

	tab_separated = query_targets.groupby(["#query", "Orthogroup"])[["Gene stable ID", "HGNC symbol", "Protein stable ID"][0]].apply("|".join).reset_index()
	for col in ["Gene stable ID", "HGNC symbol", "Protein stable ID"][1:]:
		subset = query_targets.groupby(["#query", "Orthogroup"])[col].apply("|".join).reset_index()
		tab_separated[col] = subset.loc[:, col]

	out = tab_separated.groupby([group_by_col])[non_query_cols[0]].apply(",".join).reset_index()
	for col in non_query_cols[1:]:
		subset = tab_separated.groupby([group_by_col])[col].apply(",".join).reset_index()
		out = subset.merge(out, how = "outer", on = group_by_col)

	return out 

	'''
	group_by_col = "#query"
	non_query_cols = [colname for colname in query_targets.columns.values if not colname == group_by_col]

	out = query_targets.groupby([group_by_col])[non_query_cols[0]].apply(",".join).reset_index()
	for col in non_query_cols[1:]:
		subset = query_targets.groupby([group_by_col])[col].apply(",".join).reset_index()
		out = subset.merge(out, how = "outer", on = group_by_col)

	return out 
	'''




def emapper_annotation(emapper, query_orthogroups, keep_all_targets = True):
	'''
	This function takes in emapper results and a pre-created dataframe wth the original querys and their respective orthogroups
	and combines them to output a list with all genes that share orthogroup with your query genes.
	'''


	match_column = "eggNOG_OGs"
	ortho_cols = [colname for colname in query_orthogroups.columns.values if colname.startswith("Orthogroup")]
	targets_with_orthogroups = pd.DataFrame()
	
	emapper = emapper[["#query", "eggNOG_OGs"]]
	emapper.dropna(inplace = True) # remove last three lines with emapper run data. The rest have "-" instead of NAs so we are not loosing anything

	# Find othogroups matches between query_orthogroups and all of the genes of our target species
	for col in ortho_cols:
		tax_level = col.replace("Orthogroup", "")
		subsetted_emapper = query_table(emapper, match_column = match_column, taxID = tax_level, data_origin = "emapper")
		targets_with_orthogroups = targets_with_orthogroups.append(subsetted_emapper)
	
	targets_with_orthogroups[match_column] = targets_with_orthogroups[match_column].str.replace("\|.*$", "", regex = True)
	
	# We have query_orthogroups with all orthorgoups of our curated list, and targets_with... with all target genes that have an rthorgoup in the same level at least
	# We merge them
	query_orthogroups = format_quer_orth(query_orthogroups, ortho_cols)
	query_targets= query_orthogroups.merge(targets_with_orthogroups, how = "right", left_on="Orthogroup", right_on=match_column)
	query_targets = query_targets.drop(columns = [match_column, "index"])
	if not keep_all_targets:
		query_targets = query_targets.dropna()

	out = format_query_targets(query_targets)

	return out





## Script I would put in a differnt file to make the final main function
# Import data
eggnog = read_eggnog('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Bilateria(33213)_members.tsv', '/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Metazoa(33208)_members.tsv')
lookup = pd.read_csv('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')
query = pd.read_csv("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", sep = "\t")
emapper = pd.read_csv("/g/arendt/Javier/Python/TF_annot_methods/Capitella_teleta/Data/Capitella_teleta_Prot_emapper_annotations.txt", skiprows=4, sep="\t")


# make |query - orthogroup| table
translated_eggnog = egg_translate(eggnog, lookup)
#translated_QC(translated_eggnog, query)
query_orthogroups = merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
# Get query orthogroup matching from target species genes
annotated_genes = emapper_annotation(emapper, query_orthogroups, keep_all_targets= False)
annotated_genes.to_csv("/g/arendt/Javier/Python/geneannotator/tests/Ortho_method_Capitella_TFs.tsv", sep = "\t")