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
                         names=["X", "Orthogroup", "N_Prots", "N_Spec", "ProtID", "SpeciesID"]
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
    prot_ortho = pd.DataFrame()

    eggnog["ProtID"] = eggnog.ProtID.str.split(',')
    for i in eggnog.index.values:
        prot = [x for x in eggnog.loc[i, "ProtID"] if x.startswith(taxID)]
        if len(prot) > 0:
            prot = ",".join(prot)            
            prot = pd.DataFrame(data = {"Orthogroup":[eggnog.loc[i, "Orthogroup"]], "ProtID": [prot]})
            prot_ortho = prot_ortho.append(prot, ignore_index = True)
    
    prot_ortho.drop_duplicates(inplace = True)

    if remove_taxid:
        prot_ortho.ProtID = prot_ortho.ProtID.str.replace(taxID + ".", "", regex=False)

    if explode:
        prot_ortho.ProtID = prot_ortho.ProtID.str.split(",")
        prot_ortho = prot_ortho.explode("ProtID")

    return prot_ortho


def egg_translate(eggnog, lookup, taxID = "9606", merge = True):
    """
    Takes in one or several Eggnog database raw datasets , finds each of the human proteins in it and creates a translates each one of them to Ensembl Gen IDs and HGNC symbols.
    the output is a dataframe with all human genIDs in Eggnog, their HGNC conversion and their orthogroup for all eggnog tax levels chosen
    This process takes around 20 minutes per eggnog dataset, so it is recommended to output this information in a desired folder, to avoid repeating this process
    ...

    Attributes
    ----------
    eggnog : array
        Array containing one or more eggnog datasets in pandas dataframe format.
    lookup : pandas.dataframe
        DataFrame containing three columns: "HGNC symbol", "Gene stable ID", "Protein stable ID". Each column contain strings with ID conversions from HGNC to Ensembl GenID to Ensembl protein ID
    merge : boolean
        if true, when having more than one eggnog database, it will merge all dataframes by protein ID (keeping all data from all datasets). Otherwise it outputs an array with each separated ttranslated table.
    """

    
    # Data preparation
    if not isinstance(eggnog, list):
        eggnog = [eggnog]

    for df in list(range(len(eggnog))):
        eggnog[df] = eggnog[df][eggnog[df].SpeciesID.str.contains(taxID)] # only rows with human prots stay.
        eggnog[df].drop(columns=["SpeciesID", "X", "N_Prots", "N_Spec"], inplace=True)

    lookup.dropna(subset=['Protein stable ID', "Gene stable ID"], inplace=True)


    ## Make translated table(s)
    dfs = []
    for egg in eggnog:
        egg_prots = eggnog_orthoprot_table(eggnog = egg, taxID = taxID)
        egg_prots = egg_prots.merge(lookup, how = "left", right_on = "Protein stable ID", left_on="ProtID")
        dfs.append(egg_prots)
    
    if len(dfs) > 1:
        #Iterative merge
    else: dfs = dfs[0]

    return dfs


eggnog = read_eggnog('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Bilateria(33213)_members.tsv', '/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Eggnog_Metazoa(33208)_members.tsv')
lookup = pd.read_csv('/g/arendt/Javier/Python/Human_TF_Orthogroups/TF_Data/Biomart_Lookup_Prot-HGNC-Gen_Translate_Updated.txt', sep='\t')

print(egg_translate(eggnog, lookup))