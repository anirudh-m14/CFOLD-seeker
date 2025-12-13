# CFOLD-seeker
This repo contains the CFOLD-seeker package


Dependencies:
  * os (included with python >3.7)
  * uuid (included with python >3.7)
  * time (included with python >3.7)
  * requests
  pathlib (included with python >3.7)
  concurrent.futures (included with python >3.7)
  Bio
  xml.etree.ElementTree (included with python >3.7)
  requests (included with python >3.7)
  pandas
  numpy
  ast (included with python >3.7)
  sys (included with python >3.7)

Installation and use of CFOLD-seeker:
# before using this package, be sure that those libraries are installed
  pip install requests pandas numpy biopython

Download the CFOLD-seeker package.
Inzip the package and put it in your prefered path.
Don't remove any from the existing files from this packagee.
Change in the command line to the path where CFOLD-seeker is located.
write the next line in your command line to start using CFOLD-seeker:
python CFoldSeekerCMD.py -h
This should give the next code:
"""
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
  
