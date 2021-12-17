import goterms
import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Outputs all genes/proteins in an emapper output that have a specific GO:Term in the 'GOs' column"
    )
    parser.add_argument(
        "-em",
        "--emapper",
        type=str,
        metavar="",
        required=True,
        help="Path to output file from emapper, product of running emapper on a set of proteins or a full proteome",
    )
    parser.add_argument(
        "-g",
        "--goterm",
        type=str,
        metavar="",
        required=True,
        help="GO:Term you want to search in emapper. For example 'GO:0003700' for transcription factors",
    )
    parser.add_argument(
        "-x",
        "--extra_columns",
        action="append",
        type=str,
        metavar="",
        help="One or several names of columns in emapper that you want to keep in the final output. Should have a -x flag per extra column. Default, only '#query', the names of the genes",
    )
    parser.add_argument(
        "--keep_all_columns",
        action="store_true",
        help="Default behaviour: keep only the '#query' column plus the columns specified in --extra_columns flag",
    )
    args = parser.parse_args()
    if args.keep_all_columns:
        keep_all_columns = True
    else:
        keep_all_columns = False

    if args.extra_columns is None:
        extra = False
    else:
        extra = args.extra_columns

    result = goterms.GOTerms_annotation(
        emapper=args.emapper,
        extra_columns=extra,
        GOterm=args.goterm,
        keep_all_columns=keep_all_columns,
    ).to_csv(sep="\t", index=False)
    print(result)
