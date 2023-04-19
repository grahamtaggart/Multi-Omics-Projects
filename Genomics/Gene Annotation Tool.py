from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests
import time
from lxml import etree

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

    def go_analysis(self):
        pro_seq = [str(orf[3]) for orf in self.orfs]

        if self.testing:
            pro_seq = pro_seq[:10]

        for sequence in pro_seq:
            blast_url = "https://www.uniprot.org/blast/uniprot/"
            blast_data = {
                "sequence": sequence,
                "database": "uniprotkb",
                "stype": "protein",
                "format": "out",
                "email": "graham.n.taggart.dut@gmail.com"
            }

            blast_response = requests.post(blast_url, data=blast_data)

            if blast_response.status_code == 200:
                try:
                    # Extract the job ID from the BLAST response
                    job_id = blast_response.text.strip()

                    # Check the status of the BLAST job and wait until it's finished
                    status_url = f"https://www.uniprot.org/blast/{job_id}.status"
                    while True:
                        status_response = requests.get(status_url)
                        if status_response.text.strip() == "RUNNING":
                            time.sleep(5)
                        elif status_response.text.strip() == "FINISHED":
                            break
                        else:
                            print(f"Unexpected status: {status_response.text.strip()}")
                            break

                    # Retrieve the BLAST results
                    result_url = f"https://www.uniprot.org/blast/{job_id}.out"
                    result_response = requests.get(result_url)

                    if result_response.status_code == 200:
                        try:
                            # Parse BLAST XML output to get UniProt IDs using lxml
                            root = etree.fromstring(result_response.text)
                            hits = root.xpath(".//hit")
                            uniprot_ids = [hit.find("id").text for hit in hits]

                            for uniprot_id in uniprot_ids:
                                annotation_url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
                                annotation_response = requests.get(annotation_url)

                                if annotation_response.status_code == 200:
                                    print(f"Annotations retrieved for UniProt ID: {uniprot_id}, Sequence: {sequence}")
                                else:
                                    print(
                                        f"Failed to retrieve annotations for UniProt ID: {uniprot_id}, Sequence: {sequence}")
                        except etree.XMLSyntaxError:
                            print(f"Failed to parse XML results for Sequence: {sequence}")
                            print(f"XML content: {result_response.text}")
                    else:
                        print(f"Failed to retrieve BLAST results for Sequence: {sequence}")
                        print(f"Status code: {result_response.status_code}")
                        print(f"Response content: {result_response.text}")
                except ValueError:
                    print(f"Failed to extract job ID for Sequence: {sequence}")
                    print(f"Response content: {blast_response.text}")
            else:
                print(f"BLAST job failed for Sequence: {sequence}")
                print(f"Status code: {blast_response.status_code}")
                print(f"Response content: {blast_response.text}")

            # Time delay because I'm afraid of being banned by Uniprot
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

