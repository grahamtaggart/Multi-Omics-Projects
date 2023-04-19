from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class FastaParser():
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_fasta_file(self):
        with open(self.file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record

# The GeneInfo class will take the sequence and then run a couple functions through it, calculating gc content, checking
# for open reading frames manually, then using a gene ontology API to see what's there

class GeneInfo():
    def __init__(self, seq_record: SeqRecord):
        self.seq = seq_record.seq
        self.seq_record = seq_record

    def gc_content(self) -> float:
        g = c = a = t = n = 0
        for base in self.seq:
            if base == "G":
                g += 1
            elif base == "C":
                c += 1
            elif base == "A":
                a += 1
            elif base == "T":
                t += 1
            elif base == "N":
                n += 1
            else:
                pass
        total = g + c + a + t + n
        gc_calc = ((g + c) / total) * 100
        print(f"The GC Content is: {gc_calc:.2f}%")
        return gc_calc

    def find_orfs_with_trans(self, trans_table: int, min_protein_length: int):
        answer = []
        seq_len = len(self.seq)
        for strand, nuc in [(+1, self.seq), (-1, self.seq.reverse_complement())]:
            for frame in range(3):
                trans = nuc[frame:].translate(trans_table)
                trans_len = len(trans)
                aa_start = 0
                aa_end = 0
                while aa_start < trans_len:
                    aa_end = trans.find("*", aa_start)
                    if aa_end == -1:
                        aa_end = trans_len
                    if aa_end - aa_start >= min_protein_length:
                        if strand == 1:
                            start = frame + aa_start * 3
                            end = min(seq_len, frame + aa_end * 3 + 3)
                        else:
                            start = seq_len - frame - aa_end * 3 - 3
                            end = seq_len - frame - aa_start * 3
                        answer.append((start, end, strand, trans[aa_start:aa_end]))
                    aa_start = aa_end + 1
        answer.sort()
        return answer

    def print_orfs(self, orf_list):
        for start, end, strand, pro in orf_list:
            print(
                "%s...%s - length %i, strand %i, %i:%i"
                % (pro[:30], pro[-3:], len(pro), strand, start, end)
            )

    def go_analysis(self):

# Usage example:
fasta_file = "C:\\Users\\graha\\Downloads\\ecoli.fasta"
sequence = FastaParser(fasta_file)

for i, record in enumerate(sequence.parse_fasta_file()):
    if i >= 10:  # Stop after processing the first 10 rows
        break

    print(f"Record: {record}")
    print(f"Description: {record.description}\t")
    print(f"Sequence length: {len(record.seq)}\t")

    gene_info = GeneInfo(record)
    gene_info.gc_content()

    trans_table = 11
    min_protein_length = 100
    orf_list = gene_info.find_orfs_with_trans(trans_table, min_protein_length)
    gene_info.print_orfs(orf_list)
    print()

