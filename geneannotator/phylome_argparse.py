import phylome
import argparse
import pandas as pd
import numpy as np
import os


def main(args):
    # Default value for suffix
    if args.suffix == None:
        args.suffix = "_annotated_orthology"

    ### Script
    # IF HGNC, run HGNC and exit
    if args.HGNC_method:
        print("Running HGNC method")
        annotated_tables = phylome.annotate_orthology_HGNC_method(
            args.query, args.ortho_tables
        )
        phylome.save_annotated(annotated_tables, args.output)
        exit()

    ## If user inputs translated tables, bypass making lookup and translating, just annotate with query
    elif isinstance(args.translated_orthotables, str):
        if os.path.isfile(args.translated_orthotables):

            translated_orthologies = list(
                pd.read_csv(args.translated_orthotables, sep="\t")
            )
            annotated_tables = phylome.find_query_orthologs(
                args.query, translated_orthologies
            )

        else:

            annotated_tables = phylome.find_query_orthologs(
                args.query, args.translated_orthotables
            )  # If a directory, the function will take care of it

        phylome.save_annotated(annotated_tables, args.output, args.suffix)

    else:  # --translated_orthotables not inputted

        # If user does not input translated tables but inputs a lookup table
        if isinstance(args.lookup, str):

            lookup = pd.read_csv(args.lookup, sep="\t", keep_default_na=False)
            translated_orthologies = phylome.translate_orthologies(
                args.ortho_tables, lookup
            )

            # If user inputs a lookup table and a --output_translated_orthologies
            if args.output_translated_orthotables:

                if args.suffix == "_annotated_orthology":
                    args.suffix = "_translated"

                phylome.save_annotated(translated_orthologies, args.output, args.suffix)
                exit(
                    "Translated orthology tables created in "
                    + args.output
                    + ". Now proceed to run the function with the --translated_orthotables flag and remove the --output_output_translated_orthotables"
                )

            # If doesn't input the translated tables, annotate genes
            annotated_tables = phylome.find_query_orthologs(
                args.query, translated_orthologies
            )
            phylome.save_annotated(annotated_tables, args.output, args.suffix)

        # If user does not input a lookup or translated_tables but the minimum is there (orthotables + query) run everything
        elif args.lookup is None:

            lookup = phylome.make_lookup(args.ortho_tables)

            # IF user inputs a path to save the newly made lookup
            if args.output_lookup:

                lookup.to_csv(args.output + "lookup.tsv", index=False, sep="\t")
                exit(
                    "Lookup table created in "
                    + args.output
                    + ". Now proceed to un the function with the --lookup flag and remove the --output_lookup"
                )

            translated_orthologies = phylome.translate_orthologies(
                args.ortho_tables, lookup
            )

            # If user does not input lookup or translated_tables or path to output lookup but they put --output_translated_orthotables..,
            if args.output_translated_orthotables:

                phylome.save_annotated(translated_orthologies, args.output, args.suffix)
                exit(
                    "Translated orthology tables created in "
                    + args.output_translated_orthotables
                    + ". Now proceed to run the function with the --translated_orthotables flag and remove the --output_output_translated_orthotables"
                )

            annotated_tables = phylome.find_query_orthologs(
                args.query, translated_orthologies
            )
            phylome.save_annotated(annotated_tables, args.output, args.suffix)

        else:
            exit(
                "If you want to input a lookup, please write a correct path to lookup file"
            )


if __name__ == "__main__":
    ##### Arguments parser #####
    parser = argparse.ArgumentParser(
        description="Subsets phylome orthology tables to include only proteins that have as orthologs genes/proteins from your human query. Then adds somecolumns to excplicit which genes/proteins are the orthologs"
    )

    # Full pipeline
    parser.add_argument(
        "-ot",
        "--ortho_tables",
        type=str,
        metavar="DIR",
        help="Path to folder containing all the unadulterated phylome orthology files you want to annotate or a path to one of those files.",
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        metavar="FILE",
        help="Path to file with curated list of human genes representing you module/family of interest. Genes should be in ENSEMBL_GENE_ID format with default options or HGNC (genecards) format with --HGNC flag",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="",
        required=True,
        help="Path to folder where you want to save outputs. Can be saving a lookup, translated tables or annotated tables.",
    )

    # Shortcuts
    parser.add_argument(
        "-l",
        "--lookup",
        type=str,
        metavar="",
        required=False,
        help="Optional. Path to lookup table generated by this very function to shortcut the creation of a lookup table. Should contain three columns with UniProtKB, ENEMBL_ID and HGNC_symbol conversions. Please note that the lookup table is done for the input files, so make sure you are using the lookup made with same or more files than your current input",
    )
    parser.add_argument(
        "-tr",
        "--translated-orthotables",
        type=str,
        metavar="",
        required=False,
        help="Optional. Path to translated Orthology tables. This flag will shortcut the creation of a lookup table and the translation of the orthology files.",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        type=str,
        metavar="",
        required=False,
        help="Optional. name of output files will be <taxID><suffix>.tsv . Default '_annotated_orthology' for annotated files, '_translated' for translated files. Use it if you are outputting the annotated files (end of pipeline) or --output_translated_orthotables",
    )

    # Output extra files
    parser.add_argument(
        "--output-lookup",
        action="store_true",
        help="Optional. Path to file where you want to store the newly created lookup. The generated lookup will be done specifically for the inutted ortholohy tables, so make that when you use this lookup table you do so with the same (or more) orthology tables than the ones used to make this lookup",
    )
    parser.add_argument(
        "--output-translated-orthotables",
        action="store_true",
        help="Optional. Path to folder where you want to store the newly created translated orthology tables. This translation consists of adding an extra column with all human orthologs in ENSEMBL gen ID format (instead of UNIPROT)",
    )

    # HGNC method
    parser.add_argument(
        "--HGNC-method",
        action="store_true",
        help="Use HGNC to do the matching. This avoids making a lookup and translating the orthology tables",
    )

    args = parser.parse_args()
    main(args)
