import sys
import os
import Cfoldseeker1
import Sort_csv
import cluster

HELP_TEXT = """
Usage:
    python CFoldSeekerCMD.py <input_path> [options]

Description:
    This program runs a 3-step pipeline:
        1. Cfoldseeker1.py uses Foldseek to find homologs of the input pdb files 
            and find the coordinates, assembly_IDs and scaffold_IDs of those homologous genes using Entrez-Direct.
        2. Sort_csv.py sorts all the genes and gives an intermediate csv file as output.
        3. cluster.py runs through the requirements (distance between genes, length of the whole cluster and the 
            minimum number of genes) to be a cluster and gives the final resulting clusters.


Options:
    -h, --help              Show this help message.
    --output=PATH           Give the output file path. If not defined by the user, the output will be generated in the 
                                working directory.
    --mapping=PATH          Give the mapping file path which is used to run the ID mapping from Uniprot to NCBI locally.
                                If not defined by the user, the mapping file is expected to be in the working directory.
                                The basic mapping file is part of the CFoldSeeker package.
    --evalue=X              Give a numerical value for the maximum evalue used for Foldseek. If not defined by the user
                                , 1e-5 is taken as the default value.
    --distance=X            Give a numerical value for the maximum distance between two genes of the cluster. If not 
                                defined by the user, 5000 is taken as the default value.
    --length=X              Give a numerical value for the maximum length of the whole cluster. If not defined by the 
                                user, 10000 is taken as the default value.
    --numberOfGenes=X              Give a numerical value for the minimum number of genes of the cluster. If not defined by 
                                the user, 3 is taken as the default value.
"""


def interface(argv):
    if len(argv) < 2:
        return None, {}
    input = argv[1]
    input_path = sys.argv[0]
    args = {}

    if not os.path.isfile(input_path):
        print(f"Error: File '{input_path}' does not exist.")
        sys.exit(1)
    mapping_file = sys.argv[1]
    if not os.path.isfile(mapping_file):
        print(f"Error: File '{mapping_file}' does not exist.")
        sys.exit(1)


def parse_options(args):
    opts = {
        "mapping": "idmapping_geneid.tsv",
        "output": "",
        "evalue": "1e-5",
        "distance": "5000",
        "length": "20000",
        "numberOfGenes": "3"
    }

    for arg in args:
        if arg.startswith("--"):
            if "=" in arg:
                key, val = arg[2:].split("=", 1)
                if key in opts:
                    opts[key] = val
                else:
                    print(f"Warning: unknown option '{key}'")
            else:
                print(f"Warning: ignoring malformed option '{arg}'")

    return opts


def main():
    args = sys.argv
    if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help"):
        print(HELP_TEXT)
        return
    elif len(sys.argv) < 1:
        print("Error: Missing input file.\n")
        print(HELP_TEXT)
        return
    elif len(sys.argv) >= 2:
        input = sys.argv[1]
        print(input)
        options = parse_options(sys.argv[2:])
        out1 = Cfoldseeker1.run(input,
                                mapping=options["mapping"],
                                evalue=options["evalue"])
        print("Coordinates retrieved, starting to sort")
        out2 = Sort_csv.run()
        out3 = cluster.run(output=options["output"],
                           distance=options["distance"],
                           length=options["length"],
                           numberOfGenes=options["numberOfGenes"])

if __name__ == "__main__":
    main()