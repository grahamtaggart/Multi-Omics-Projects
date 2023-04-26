from Bio import SeqIO
from colorama import Fore, Style


# By using this parse_fasta_file func I can use as many different fasta files as I want without refactoring SeqIO.parse everytime

def parse_fasta_file(file_path):
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            yield record


for record in parse_fasta_file(r"C:\Users\graha\Downloads\ecoli.fasta"):
    sequence = record

# Probably some faster algo can do this without major performance issues, but I don't know about it yet.

def colors_to_letters(sequence: str) -> str:
    for base in sequence:
        if base == "A":
            print(Fore.LIGHTRED_EX + "A" + Style.RESET_ALL, end="")
        elif base == "T":
            print(Fore.LIGHTYELLOW_EX + "T" + Style.RESET_ALL, end="")
        elif base == "G":
            print(Fore.LIGHTGREEN_EX + "G" + Style.RESET_ALL, end="")
        elif base == "C":
            print(Fore.LIGHTBLUE_EX + "C" + Style.RESET_ALL, end="")
        elif base == "N":
            print(Fore.LIGHTCYAN_EX + "N" + Style.RESET_ALL, end="")
        else:
            print("")

def main():
    colored_seq = colors_to_letters(sequence)
    print(colored_seq)

if __name__ == "__main__":
    main()
