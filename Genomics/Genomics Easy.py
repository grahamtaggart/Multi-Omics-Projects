from Bio import SeqIO
from colorama import Fore, Style
import argparse

def colors_to_letters(seq: str):
    for base in seq:
        if base == "A":
            print(Fore.LIGHTRED_EX + "A" + Style.RESET_ALL, end="")
        elif base == "T":
            print(Fore.LIGHTYELLOW_EX + "T" + Style.RESET_ALL, end="")
        elif base == "G":
            print(Fore.LIGHTGREEN_EX + "G" + Style.RESET_ALL, end="")
        elif base == "C":
            print(Fore.LIGHTBLUE_EX + "C" + Style.RESET_ALL, end="")

seq = "ACGTCGTCGTCAACACACACACACACACAGTGTGGTGTGTGTGTGCCACACACAGTGGTG"

colors_to_letters(seq)







