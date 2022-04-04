import sys
from os.path import exists

STOP = ['TAA', 'TAG', 'TGA']
START = 'ATG'

CODON_AA_MAPPING = \
    {"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
      "TGT": "C", "TGC": "C",
      "GAT": "D", "GAC": "D",
      "GAA": "E", "GAG": "E",
      "TTT": "F", "TTC": "F",
      "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
      "CAT": "H", "CAC": "H",
      "ATA": "I", "ATT": "I", "ATC": "I",
      "AAA": "K", "AAG": "K",
      "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
      "ATG": "M", "AAT": "N", "AAC": "N",
      "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
      "CAA": "Q", "CAG": "Q",
      "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
      "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
      "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
      "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
      "TGG": "W",
      "TAT": "Y", "TAC": "Y",
      "TAA": "_", "TAG": "_", "TGA": "_"}

COMMAND_HELP = "python run.py [FASTA_FILE] [CSV_PROTEIN_FILE]"


def main():
    validate_arguments()
    calculate_and_print_stats()


def validate_arguments():
    if len(sys.argv) != 3:
        raise ValueError("Incorrect number of arguments!\n%s" % COMMAND_HELP)

    if not exists(sys.argv[1]) or not exists(sys.argv[2]):
        raise ValueError("Invalid file(s), please verify!\n%s" % COMMAND_HELP)


def gc_content(dna_seq):
    if len(dna_seq) == 0:
        return 0

    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc":
            gc_count += 1
    return gc_count / len(dna_seq)


def update_frequency(seq, current_dic):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    for s in seq.upper():
        if s in current_dic:
            current_dic[s] += 1


def is_eof(file_line):
    return not file_line


def print_stats_results(sequence_len, nucleotides_stats, gc_content, codons_stats):
    print("1. Sequence length => %s" % sequence_len)
    print("2. Frequencies => 'A': %s, 'T': %s, 'C': %s, 'G': %s" % (
        nucleotides_stats['A'], nucleotides_stats['T'], nucleotides_stats['C'], nucleotides_stats['G']))
    print("3. GC content => %s" % gc_content)
    print("4. Number of Start codons (AUG) =>  %s" % codons_stats['AUG'])
    print("5. Number of Stop codons UAA => %s, UAG: %s, UGA: %s" % (
        codons_stats['UAA'], codons_stats['UAG'], codons_stats['UGA']))

    sorted_codons_stats = sorted(codons_stats.items(), key=lambda x: x[1])
    print("6. Most (%s) / Less (%s) frequent codons" % (
        sorted_codons_stats[len(sorted_codons_stats) - 1][0], sorted_codons_stats[0][0]))


def complement(codon):
    return codon.replace('A', 'U').replace('T', 'A').replace('C', '-').replace('G', 'C').replace('-', 'G')


def update_codons(seq_line, codons_stats):
    for idx_codon in range(0, len(seq_line), 3):
        codon = complement(seq_line[idx_codon: idx_codon + 3])
        if codon in codons_stats:
            codons_stats[codon] += 1
        elif len(codon) == 3:
            codons_stats[codon] = 1
        else:
            return seq_line[idx_codon: len(seq_line)]  # returns the incomplete codon

    return ""


def calculate_and_print_stats():
    fasta_file = open(sys.argv[1], "r", 1024)
    full_dna_seq = ""
    try:
        sequence_len = 0
        nucleotides_stats = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        codons_stats = {}
        remaining = ""
        fasta_file.readline()  # skips the first line
        while True:
            seq_line = fasta_file.readline()
            if is_eof(seq_line):
                break
            seq_line = seq_line.replace('\n', '')
            full_dna_seq += seq_line
            sequence_len += len(seq_line)
            update_frequency(seq_line, nucleotides_stats)
            remaining = update_codons(remaining + seq_line, codons_stats)
    finally:
        fasta_file.close()

    gc_count = gc_content(full_dna_seq)
    print_stats_results(sequence_len, nucleotides_stats, gc_count, codons_stats)
    find_orfs_and_persist_data(full_dna_seq)
    do_overlap()


def find_orfs_and_persist_data(full_dna_seq):
    min_nucleotides_len = 30
    orf_coordinates = get_orfs_and_coordinates(full_dna_seq, min_nucleotides_len)
    save_orf_coordinates(orf_coordinates)
    save_protein_seqs(orf_coordinates)


def save_orf_coordinates(orf_coordinates):
    file = open("orf_coordinates.txt", "w")
    try:
        for orf_coordinate in orf_coordinates:
            file.write(str(orf_coordinate["start"]) + "," + str(orf_coordinate["stop"]) + "," + orf_coordinate["orf"] + "\n")
    finally:
        file.close()


def save_protein_seqs(orf_coordinates):
    file = open("all_potential_proteins.txt", "w")
    try:
        for orf_coordinate in orf_coordinates:
            arr_possible_proteins = all_proteins_rf(orf_coordinate["aa"])
            for possible_protein in arr_possible_proteins:
                file.write(possible_protein + "\n")
    finally:
        file.close()


def codon_to_aa(codon):
    if codon in CODON_AA_MAPPING:
        return CODON_AA_MAPPING[codon]
    else:
        return None


# Computes all possible proteins in an aminoacid sequence
def all_proteins_rf(aa_seq):
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def get_index_start_codon(sequence, starting_point=0):
    for index in range(starting_point, len(sequence), 3):
        if sequence[index: index + 3] == START:
            return index

    return -1


def get_orfs_and_coordinates(sequence, min_nucleotides_len):
    orf_coordinates_data = []

    for index_frame in range(0, 3):   # emulates the positive reading-frames
        start_codon_index = get_index_start_codon(sequence, index_frame)
        frame_value = index_frame + 1

        if start_codon_index == -1:
            continue     # No start codon found into the sequence

        current_orf = START
        aa_seq = codon_to_aa(START)
        codon_index = start_codon_index + 3
        while codon_index < len(sequence):
            codon = sequence[codon_index: codon_index + 3]

            if len(codon) != 3:     # if reaches the end with a number of nt lower than 3, breaks out the loop
                break

            aa_seq += codon_to_aa(codon)
            current_orf += codon

            if codon in STOP:
                nucleotides_count = len(current_orf)
                if nucleotides_count >= min_nucleotides_len:
                    end_stop_index = codon_index + 3    # adds 3 to account for the full nt triplet + the '0' index starting
                    orf_coordinates_data.append({"start": start_codon_index + 1, "stop": end_stop_index, "orf": current_orf, "aa": aa_seq, "frame": frame_value, "nt_len": len(current_orf), "aa_len": len(aa_seq)})

                start_codon_index = get_index_start_codon(sequence, codon_index + 3)
                if start_codon_index == -1:  # no further ORFs into the remaining sequence
                    break

                codon_index = start_codon_index + 3
                current_orf = START
                aa_seq = codon_to_aa(START)
                continue

            codon_index += 3

    return orf_coordinates_data


def do_overlap():
    pass


if __name__ == "__main__":
    main()
