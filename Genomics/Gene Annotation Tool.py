from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import time
import requests

class FastaParser():
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_fasta_file(self):
        with open(self.file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record

# The GeneInfo class will take the sequence and then run a couple functions through it, calculating gc content, checking
# for open reading frames manually, then using Uniprot's APIs to check homology, then grabbing annotation data if
# homology is found.

class GeneInfo():
    def __init__(self, seq_record: SeqRecord, testing=False):
        self.seq = seq_record.seq
        self.seq_record = seq_record
        self.orfs = []
        self.testing = testing

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
        self.orfs = answer
        return answer

    def print_orfs(self, orf_list):
        for start, end, strand, pro in orf_list:
            print(
                "%s...%s - length %i, strand %i, %i:%i"
                % (pro[:30], pro[-3:], len(pro), strand, start, end)
            )

    def submit_blast_job(self, sequence):
        print("Submitting BLAST job for sequence:", sequence)
        try:
            # Using a shorter database "pdb" for testing
            result_handle = NCBIWWW.qblast("blastp", "pdb", sequence)
            print("BLAST job submitted.")
            print("Result handle type:", type(result_handle))
            print("Result handle content:", result_handle.read()[:100])  # Print the first 100 characters of the result
            result_handle.seek(0)  # Reset the result handle's position
        except Exception as e:
            print("Error occurred while submitting BLAST job:")
            print(str(e))
            return None
        return result_handle

    def parse_blast_results(self, result_handle):
        print("Parsing BLAST results.")
        try:
            blast_record = NCBIXML.read(result_handle)
            print("BLAST results parsed.")
        except Exception as e:
            print("Error occurred while parsing BLAST results:")
            print(str(e))
            return None
        return blast_record

    def get_uniprot_ids_from_blast_record(self, blast_record):
        print("Extracting UniProt IDs from BLAST record.")
        try:
            uniprot_ids = [hit.accession.split("|")[1] for hit in blast_record.alignments]
            print(f"UniProt IDs extracted: {uniprot_ids}")
        except Exception as e:
            print("Error occurred while extracting UniProt IDs from BLAST record:")
            print(str(e))
            return []
        return uniprot_ids

    def get_uniprot_annotation(self, uniprot_id):
        print(f"Retrieving UniProt annotation for ID: {uniprot_id}")
        annotation_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
        annotation_response = requests.get(annotation_url)
        print(f"UniProt annotation retrieved. Status code: {annotation_response.status_code}")
        return annotation_response

    def go_analysis(self):
        pro_seq = [str(orf[3]) for orf in self.orfs]

        if self.testing:
            pro_seq = pro_seq[:10]

        for sequence in pro_seq:
            try:
                result_handle = self.submit_blast_job(sequence)
                if result_handle is None:
                    print(f"Failed to submit BLAST job for sequence: {sequence}")
                    continue

                blast_record = self.parse_blast_results(result_handle)
                if blast_record is None:
                    print(f"Failed to parse BLAST results for sequence: {sequence}")
                    continue

                uniprot_ids = self.get_uniprot_ids_from_blast_record(blast_record)

                for uniprot_id in uniprot_ids:
                    annotation_response = self.get_uniprot_annotation(uniprot_id)

                    if annotation_response.status_code == 200:
                        print(f"Annotations retrieved for UniProt ID: {uniprot_id}, Sequence: {sequence}")
                    else:
                        print(
                            f"Failed to retrieve annotations for UniProt ID: {uniprot_id}, Sequence: {sequence}")

            except Exception as e:
                print(f"An error occurred while processing the sequence: {sequence}")
                print(str(e))

            # Time delay to avoid overloading the server
            time.sleep(5)
# Main use:

def main():
    fasta_file = "C:\\Users\\graha\\Downloads\\ecoli.fasta"
    sequence = FastaParser(fasta_file)

    for i, record in enumerate(sequence.parse_fasta_file()):
        if i >= 10:
            break

        print(f"Record: {record}")
        print(f"Description: {record.description}\t")
        print(f"Sequence length: {len(record.seq)}\t")

        gene_info = GeneInfo(record, testing=True)
        gene_info.gc_content()

        trans_table = 11
        min_protein_length = 100
        orf_list = gene_info.find_orfs_with_trans(trans_table, min_protein_length)
        gene_info.print_orfs(orf_list)
        gene_info.go_analysis()
        print()

if __name__ == '__main__':
    main()

