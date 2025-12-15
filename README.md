# CFOLD-seeker
This repo contains the CFOLD-seeker package.  
CFOLD-seeker is a tool to find homologous clusters starting from a known set of pdb or cif files which form a cluster.

### Dependencies:
  * requests
  * biopython
  * pandas
  * numpy
### Dependencies included in python >3.7:
  * os
  * uuid
  * time
  * pathlib
  * concurrent.futures
  * xml.etree.ElementTree
  * ast
  * sys

## Installation and use of CFOLD-seeker:
Before using this package, be sure to install these libraries:  
  `pip install requests pandas numpy biopython`

1. Open the latest release on the right-end side of the github page.
2. Download the CFOLD-seeker ZIP-file and the idmapping_geneid.tsv from the latest release.
3. Unzip the package to your preferred path and add the idmapping_geneid.tsv to this path.
4. Change directory to the path where CFOLD-seeker is located.
5. Use the next command to start using CFOLD-seeker:  
`python CFoldSeekerCMD.py -h`
 
This gives the next output:  

Usage:     
&emsp; python CFoldSeekerCMD.py <input_path> [options]

Description:  
    &emsp; This program runs a 3-step pipeline:  
    
1. Cfoldseeker1.py uses Foldseek to find homologs of the input pdb or cif files and find the coordinates, assembly_IDs and scaffold_IDs of those homologous genes using Entrez-Direct.  
2. Sort_csv.py sorts all the genes and gives an intermediate csv file as output.  
3. cluster.py runs through the requirements (distance between genes, length of the whole cluster and the minimum number of genes) to be a cluster and gives the final resulting clusters.

Input_path:  
&emsp; The path where the pdb or cif files of the cluster are located 

Options:  
    &emsp; -h, --help  
    &emsp; &emsp;         Show this help message.  
    &emsp; --output=PATH  
    &emsp; &emsp;     The output file path. If not defined by the user, the output will be generated in the working directory.  
    &emsp; --mapping=PATH  
    &emsp;  &emsp;   The mapping file path which is used to run the ID mapping from Uniprot to NCBI locally. If not defined by the user, the mapping file is expected to be in the working directory. The basic mapping file is part of the CFoldSeeker package.  
    &emsp; --evalue=X  
    &emsp; &emsp; A numerical value for the maximum evalue used for Foldseek. If not defined by the user, 1e-5 is taken as the default value.  
    &emsp; --distance=X  
    &emsp; &emsp; A numerical value for the maximum distance between two genes of the cluster. If not defined by the user, 5000 is taken as the default value.  
    &emsp; --length=X  
    &emsp; &emsp; A numerical value for the maximum length of the whole cluster. If not defined by the user, 10000 is taken as the default value.  
    &emsp; --numberOfGenes=X  
    &emsp; &emsp; A numerical value for the minimum number of genes of the cluster. If not defined by the user, 3 is taken as the default value.  
