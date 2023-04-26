# Comparative genomics - Create a program that compares the genomes of different species to identify conserved regions,
# unique genes, and potential functional elements.

#function takes in many different sequences, adds commas and other delimeters so that they are sorted, maybe adds them to a list
#from that list run a loop that iterates them through NCBI APIs to determine 1. conserved regions, 2. unique genes,
#3. functional elements. This program is actually easier than the gene annotation tool because it can all be done with
#one or two NCBI API calls I'd recon. Need a function that prints results, maybe to make it more complicated or
#informative use biopython to check differences, or try to visualize the results afterwards. Plan is
#1. write sequence function
#2. All the NCBI API functions (maybe a class to simplify)
#3. Biopython visualizer? Either that or some other data visualizer so I can better understand differences in genomes.#

from Bio import SeqIO
from Bio import Entrez
from Bio import AlignIO
import requests
import json

# Data Retrieval Step

def entrez_call():
    Entrez.email = "graham.n.taggart.dut@gmail.com"
    handle = Entrez.einfo()
    result = handle.read()
    print(result.keys())


# Data Pre-processing Step
