import pandas as pd
import numpy as np
import os
import urllib.parse
import urllib.request
import re
import httplib2 as http
import json
import time
from tqdm import tqdm


def translate_uniprots(orthotable, lookup):
    """
    Takes all orthology tables and translates each of the human Uniprot Orthologs into ENSEMBL IDs.
    Then it makes a new orhology table with the translations in extra columns.
    A single uniprot ID can have more than one ENSEMBL genID, so the pipeline separates ENSEMBL within a single Unirpot ID
    by "|" and ENS IDs from different Uniprots by ","

    Attributes
    ----------
    orthotable: String.
                    Path to folder with orthology tables
    lookup: pandas dataframe.
                    Two columns: "Gene stable ID", "UniProtKB Gene Name ID"
    """
    ENS_col = orthotable.columns.get_loc("ENSEMBL_ID")
    ortho_col = orthotable.columns.get_loc("orthologs")

    for row in list(range(len(orthotable.index))):
        query = orthotable.iat[row, ortho_col]  # orthologs
        query = query.replace("|", ",")
        query = query.split(",")

        for uniprotID in tqdm(query):
            uniprotID = uniprotID.replace("9606.", "")
            genIDs = lookup["ENSEMBL_ID"][
                lookup["UniProtKB"] == uniprotID
            ]  # there are several GenIDs per UniprotID. We will take all of them
            genIDs = "|".join(list(set(genIDs)))  # eliminate duplicates
            orthotable.iat[row, ENS_col] = orthotable.iat[row, ENS_col] + "," + genIDs

    return orthotable


def directory_or_file(path):
    """
    If given a directory makes full paths of each child file. If given a file just puts it in an array
    """
    if os.path.isdir(path):
        orthology_tables = os.listdir(path)
        orthology_tables = [path + file for file in orthology_tables]
    elif os.path.isfile(path):
        orthology_tables = [path]
    else:
        exit("Introduced path is neither a folder or a file")

    return orthology_tables


def human_genes_string(path):
    """
    This script outputs a string with all human genes in the orthology tables we might want to translate
    It takes all orthology tables, individualizes all Uniprot Human orthologs, puts them in a string,
    and finally it removes duplicates. This string is used later to make a lookup table using uniprot's API

    Attributes
    ----------
    path: String
            Absolut path to orthology tables
    ortho_directory: Boolean
            True if you want to run the pipeline in all phylome files in a folder (in which case the folder should contain ONLY phylome orthology files). False if path is to a single othology file.

    """
    # Is input path a file or a directory?
    orthology_tables = directory_or_file(path)

    genes = ""
    for fullpath in orthology_tables:
        orthoTable = pd.read_csv(
            fullpath, index_col=False, skiprows=[i for i in range(1, 13)], sep="\t"
        )
        orthoTable = orthoTable[orthoTable["target_species"] == "Homo sapiens"]
        orthoTable = orthoTable["orthologs"]

        orthoTable = orthoTable.str.replace("9606.", " ", regex=False)
        orthoTable = orthoTable.str.replace("|", " ", regex=False)
        orthoTable = orthoTable.str.replace(",", " ", regex=False)

        # Make string out of them
        query = orthoTable.str.cat(sep=" ")
        genes = genes + query

    # Make string with them, no duplicated genes.
    genes = genes.split(" ")
    genes = set(genes)
    genes = " ".join(genes)
    return genes


def uniprot_request(genes, from_id, to_id):
    """
    Uses Uniprot's api to make a lookup table having UniprotKB ID - EnsemblID - HGNC symbol

    ...

    Attributes
    ----------
    genes : str
        a formatted string with all Uniprot ids to translate, separated by "\t"
    from_format : str
        Uniprot ID abreviation as in https://www.uniprot.org/help/api_idmapping of the inputted genes to translate
    to_format : str
        Gene code abbreviation as in https://www.uniprot.org/help/api_idmapping to translate Unirot IDs to
    """
    url = "https://www.uniprot.org/uploadlists/"

    params = {"from": from_id, "to": to_id, "format": "tab", "query": genes}

    data = urllib.parse.urlencode(params)
    data = data.encode("utf-8")
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    return response.decode("utf-8")


def format_uniprot_output(uniprot, columns):
    uniprot = pd.DataFrame(
        [x.split("\t") for x in uniprot.split("\n")]
    )  # string to DataFrame
    uniprot.columns = columns
    end = len(uniprot.index) - 1
    uniprot = uniprot.iloc[1:end]

    return uniprot


def lost_genes(TFs, path_to_orthologies, all_lost_genes=False, TF_lost_genes=True):
    """
    this function takes in the translated orthology tables and tells you which genes were not translated if all_lost_genes=True
    If TF_lost_genes is true, then all_lost_genes must be false. In this case you will output only those lost genes taht are Tfs according to your TFs file

    Attributes
    ----------
    TFs: list.
            Object with all TFs ENSEMBL IDs
    path_to_orthologies: string
            object with path to folder with orthology tables
    """
    finalist = []
    finalUniprot = []

    for file in os.listdir(path_to_orthologies):
        fullpath = os.path.join(path_to_orthologies, file)
        orthoTable = pd.read_csv(fullpath, sep="\t")
        lostGenes = orthoTable[["GeneName_target", "orthologs"]][
            orthoTable["GenIDs"].isna()
        ]

        Uniprotlost = list(lostGenes["orthologs"])
        lostGenes = list(lostGenes["GeneName_target"])

        finalUniprot = finalUniprot + Uniprotlost
        finalist = finalist + lostGenes

    finalostable = pd.DataFrame({"Genes": finalist, "Uniprots": finalUniprot})

    if all_lost_genes:
        finalostable.drop_duplicates(subset=["Genes"], inplace=True)

    elif TF_lost_genes:
        finalist = list(set(finalist))  # unique values
        condition = np.isin(finalist, TFs)
        finalist = pd.DataFrame({"Genes": finalist})[condition]  # only TFs

        finalostable = finalostable.merge(
            finalist, how="right", left_on="Genes", right_on="Genes"
        )
        finalostable = finalostable.drop_duplicates(subset=["Genes"])
    return finalostable


def HGNC_request(gene):

    headers = {
        "Accept": "application/json",
    }

    uri = "http://rest.genenames.org"
    path = "/fetch/symbol/" + gene

    target = urllib.parse.urlparse(uri + path)
    method = "GET"
    body = ""

    h = http.Http()

    response, content = h.request(target.geturl(), method, body, headers)

    if response["status"] == "200":
        try:
            data = json.loads(content)
            return data["response"]["docs"][0]["ensembl_gene_id"]
        except Exception:
            pass

    else:
        print("Error detected: " + response["status"])
    # time.sleep(0.01) # 20 miliseconds per request


def update_Biomart(lookup, finalostable):
    """
    This function takes a the product of lost_genes() and uses it to update the Biomart lookup table

    Attributes
    ----------
    lookup: Pandas dataframe
            Biomart lookup table with threee columns ["Gene stable ID", "HGNC symbol", "UniProtKB Gene Name ID"]
    Finalostable: Pandas dataframe
            Product of lost_genes(). Two columns, one called "Genes" with HGNC symbols and one called "Uniprots" with uniprots (can be preceeded by "9606.")
    """
    ## Use final list of lost genes to update Biomart table and avoid loosing those genes
    lookup2 = lookup.drop(columns=["UniProtKB Gene Name ID"])
    lookup2.drop_duplicates(subset=["HGNC symbol"], inplace=True)
    finalostable = finalostable.merge(
        lookup2, how="left", left_on="Genes", right_on="HGNC symbol"
    )

    # Format table to fit Biomart standarst and plug it at the end of the Biomart lookup table
    finalostable.drop(columns=["Genes"], inplace=True)
    finalostable = finalostable[["Gene stable ID", "HGNC symbol", "Uniprots"]]
    finalostable["Uniprots"] = finalostable["Uniprots"].str.replace("9606.", "")
    finalostable.columns = ["Gene stable ID", "HGNC symbol", "UniProtKB Gene Name ID"]
    print("Lines that will be added to updated version of Biomart: \n")
    print(finalostable)

    lookup = lookup.append(finalostable, ignore_index=True)
    return lookup


def find_position(row, gene, column="ENSEMBL_ID", HGNC=False):
    if HGNC:
        # This regex is only needed for HGNCs. However leads to certain problems, as a HGNC in the table gives multiple matches. This is solved further in the HGNC pipeline, but we bypass it for the ENSGo one.
        position = row[column].str.split(",")
        for i in position.index.values:
            position[i] = str(list(np.isin(position[i], gene)))
            position[i] = position[i].split(
                "True"
            )  # The problem is that if you use regex to split you end up taking the comma aswell, and you cannot just use gene because many HGNCs are substrings of other, so the split gets mistaken. What we are doing is transform into a string of true false and split by True

    else:
        position = row[column].str.split(gene)

    position = position.str[0]
    position = position.str.count(
        ","
    )  # if there are 3 "," before the found ID, then the found ID is in the 4th position (python counts from 0, so no need to +1)

    return position


def query_position_table(rows, position, gene, HGNC=False):
    """
    This function takes a row from the translated orthotable that has a TF (object called rows),
    The position within the GenID string, in which our gene of interest is, and makes a dataframe containing
    for that row, which GenIDs are TFs, their HGNC tranlslation (found based on poistion in string) and the position.

    Attributes
    ----------
    rows: pandas dataframe
            containing "##Seed_(co-)orthologs", "type", "GenIDs", "Ordered_HGNCs", "GeneName_target" as columns. Only rows that contain genes that are TFs.
    position: pandas dataframe
            Contains for each of the rows in row object, in which position the TF is. This is, if GenIDs are: ENS1,ENS2,ENS3 and ENS2 is a TF, then it is in postion 2.
    gene: string.
            ENS ID of thegene we are looking for, it's just to avoid getting also the other gen IDs separated by "|"
    HGNC: boolean.
            If you are not using tranlated orthotables and just HGNCs, then True.
    """

    table = pd.DataFrame(
        {"GenID": [], "HGNC": [], "position_from_0": [], "number_of_IDs": []}
    )

    HGNC_symbols = rows["GeneName_target"].str.split(",")
    if not HGNC:
        ENS_symbols = rows["ENSEMBL_ID"].str.split(",")

    indexes = HGNC_symbols.index.values
    for i in indexes:
        position_tmp = position[i]
        HGNC_name = HGNC_symbols[i][position_tmp]

        # DISCLAIMER, this function assumees that you only have one GenID per TF. so within a srting of "|" only one geneID is correct.
        if not HGNC:
            # substitute this following chunk for ENS_name = gene to ignore all the other ENSIDs that share UnirptoID
            ENS_name = ENS_symbols[i][position_tmp]
            ENS_name = ENS_name.split("|")
            ENS_name = np.isin(ENS_name, gene)
            ENS_name = [gene if i == True else "-" for i in ENS_name]
            ENS_name = "|".join(ENS_name)
        else:
            ENS_name = ""

        number_of_IDs = len(
            HGNC_symbols[i]
        )  # total number of symbols translated or not
        new = pd.DataFrame(
            {
                "GenID": [ENS_name],
                "HGNC": [HGNC_name],
                "position_from_0": [position_tmp],
                "number_of_IDs": [number_of_IDs],
            }
        )
        table = table.append(new)

    table.index = HGNC_symbols.index
    return table
