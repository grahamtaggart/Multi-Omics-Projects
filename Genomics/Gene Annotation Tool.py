from Bio import SeqIO, SeqRecord, SeqUtils
import numpy as np
import pandas as pd
import scipy

class FastaParser:
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_fasta_file(self):
        with open(self.file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record


#file_path = "C:\\Users\\graha\\PycharmProjects\\Multi-Omics Projects\\Genomics\\Files\\spermwhalesequence.fasta"
#parser = FastaParser(file_path)

#for record in parser.parse_fasta_file():
    #sequence = record

class GeneRetrieve:
    def __init__(self, call, compare, download):
        self.call = call
        self.compare = compare
        self.download = download

    @property.setter
    def CallBlast(self, value):
        x = value + 1
        self._CallBlast = x

class GeneAnnotator:
    def __init__(self, identity, annotate, predict):
        self.identity = identity
        self.annotate = annotate
        self.predict = predict

