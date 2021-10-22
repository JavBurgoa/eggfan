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
    prot_ortho = pd.DataFrame()

    eggnog[prot_column] = eggnog.loc[:, prot_column].str.split(',')
    for i in eggnog.index.values:
        prot = [x for x in eggnog.loc[i, prot_column] if x.startswith(taxID)]
        if len(prot) > 0:
            prot = ",".join(prot)            
            prot = pd.DataFrame(data = {"Orthogroup":[eggnog.loc[i, "Orthogroup"]], prot_column: [prot]})
            prot_ortho = prot_ortho.append(prot, ignore_index = True)
    
    prot_ortho.drop_duplicates(inplace = True)

    if remove_taxid:
        prot_ortho.loc[:, prot_column] = prot_ortho.loc[:, prot_column].str.replace(taxID + ".", "", regex=False)

    if explode:
        prot_ortho.loc[:, prot_column] = prot_ortho.loc[:, prot_column].str.split(",")
        prot_ortho.loc[:, prot_column] = prot_ortho.loc[:, prot_column].explode(prot_column)

    return prot_ortho


def egg_translate(eggnog, lookup, taxID = "9606", suffixes = False):
    """
    Takes in one or several Eggnog database raw datasets , finds each of the proteins in each eggnog dataset and creates a table specifying for each protein their respective Ensembl gene ID, HGNC symbol (HGNC and Ensembl gen ID may not be "translated", Na for empty values) and orthogroup.
    ...

    Attributes
    ----------
    eggnog : list or pandas dataframe
        list containing one or more eggnog datasets in pandas dataframe format or a single pandas dataframe. Supposed to be direct output from read_eggnog.
    lookup : pandas.dataframe
        DataFrame containing three columns: "HGNC symbol", "Gene stable ID", "Protein stable ID". Each column contain strings with ID conversions from HGNC to Ensembl GenID to Ensembl protein ID
    suffixes : list
        Only should be added if inputting more than one eggnog dataset (in a list). List containing a string per eggnog dataset. When having more than one eggnog database, if there is a list, it will merge all dataframes by protein ID (keeping all data from all datasets) and adding the suffixes in the list to the orthorgoups column. If nothing specified merging will not happen
    
    Output
    ------
    - If input is a single eggnog dataset the output is a single table 
    - If input is a list with multiple eggnog datasets + suffixes argument then the output is a single table with an extra column per dataset specifying orthorgoups from the differnt datasets inputed
    """

    prot_column = "Protein stable ID"

    # Data preparation
    if not isinstance(eggnog, list):
        eggnog = [eggnog]

    for df in list(range(len(eggnog))):
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
            out = out.merge(df, on = prot_column, suffixes = suffixes)
    return out



def translated_QC(translated_eggnog, query = False, match_column = "Gene stable ID"):
    ortho_cols = [colname for colname in translated_eggnog.columns.values if colname.startswith("Orthogroup")]

    translated_eggnog = translated_eggnog.drop(columns = ortho_cols)
    lost = translated_eggnog.loc[translated_eggnog[match_column].isna(), :]
    lost_prots = lost["Protein stable ID"]
    print("Proteins thatwere in eggnog but were not translated:")
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
        print(columns_keep)
        query_with_orth = query_with_orth.loc[:, columns_keep]

    return query_with_orth







## Script I would put in a differnt file to make the final main function
eggnog = read_eggnog('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Bilateria(33213)_members.tsv', '/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Metazoa(33208)_members.tsv')
lookup = pd.read_csv('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')
query = pd.read_csv("/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/TFs_Ensembl_v_1.01.txt", sep = "\t")


translated_eggnog = egg_translate(eggnog, lookup, suffixes = ["Bilatera", "Metazoa"])
translated_QC(translated_eggnog, query)
query_orthogroups= merge_with_query(translated_eggnog, query, merge_on = "Gene stable ID", keep_conversions= True)
