from Bio import SeqIO

class FastaParser():
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_fasta_file(self):
        with open(self.file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record
fasta_file = "C:\\Users\\graha\\PycharmProjects\\Multi-Omics Projects\\Genomics\\Files\\spermwhalesequence.fasta"
sequence = FastaParser(fasta_file)

for record in sequence.parse_fasta_file():
    print(f"Record: {record}")
    print(f"Description: {record.description}\t")
    print(f"Sequence length: {len(record.seq)}\t")

# The GeneInfo class will take the sequence and then run a couple functions through it, calculating gc content, open reading frames, whether transposons are present or not,

class GeneInfo():
    def __init__(self, seq):
        self.seq = seq

    def gc_content(self):
        g = 0
        c = 0
        a = 0
        t = 0
        n = 0
        for base in sequence:
            if base == "G" or "C":
                g += 1
                c += 1
            elif base == "T" or "A":
                a += 1
                t += 1
            elif base == "N":
                n += 1
            else:
                pass
        total = g + c + a + t + n
        gc_calc = ((g + c) / total ) * 100
        print(f"The GC Content is:{gc_calc}")
    def ofr(self):
    def dna_trans(self):
    def retro_trans(self):
    def ontology(self):
