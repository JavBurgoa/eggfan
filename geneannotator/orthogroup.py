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
    Takes in one or several Eggnog database raw datasets , finds each of the human proteins in it and creates a translates each one of them to Ensembl Gen IDs and HGNC symbols.
    the output is a dataframe with all human genIDs in Eggnog, their HGNC conversion and their orthogroup for all eggnog tax levels chosen.
    3 possible inputs/outputs:
        - A single eggnog dataset / a single tranlated table
        - An list with multiple datasets + merge argument / A single translated tabe with a orthorgoup column per eggnog database
        - An list with multiple datasets + no merge argument / An list with a tranlated table per eggnog dataset introduced
    ...

    Attributes
    ----------
    eggnog : list
        list containing one or more eggnog datasets in pandas dataframe format.
    lookup : pandas.dataframe
        DataFrame containing three columns: "HGNC symbol", "Gene stable ID", "Protein stable ID". Each column contain strings with ID conversions from HGNC to Ensembl GenID to Ensembl protein ID
    suffixes : list
        list containing a string per eggnog dataset. When having more than one eggnog database, if there is a list, it will merge all dataframes by protein ID (keeping all data from all datasets) and adding the suffixes in the list to the orthorgoups column. If nothing specified merging will not happen
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
        if isinstance(suffixes, list):
            for df in dfs[1:]:
                df = df.drop(columns = ["HGNC symbol", "Gene stable ID"])
                out = out.merge(df, on = prot_column, suffixes = suffixes)
        else: out = dfs

    return out










## Script
eggnog = read_eggnog('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Bilateria(33213)_members.tsv', '/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Metazoa(33208)_members.tsv')
lookup = pd.read_csv('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')

print(egg_translate(eggnog, lookup, suffixes = ["Bilatera", "Metazoa"]))