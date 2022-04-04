# -*- coding: utf-8 -*-

def get_codon_amino_dict(filename):
  with open(filename, "r") as f:
    file_content = f.readlines()

  return dict(
      [i.replace('\n','').replace('"','').split(' ') for i in file_content]
      )
# dict_codons_to_amino = get_codon_amino_dict(path_to_folder + "CLASS_02/code_and_data/" + "/genetic_code.txt")

def start_codon_position(seq, reading_frame, k=3):
  """
  Finds the position of the start codons.
  reading_frame: can be 0, 1 or 2.
  """
  lt_init_codons = ["ATG"]
  seq = seq.upper()
  lt = []
  for i in range(reading_frame, len(seq)-k+1, k):
    if seq[i:i+k] in lt_init_codons:
      lt.append((i, seq[i:i+k]))
  if lt:
    return lt
  else:
    return False
      
def end_codon_position(seq, reading_frame, k=3):
  """
  Finds the position of the end codons.
  reading_frame: can be 0, 1 or 2.
  """
  seq = seq.upper()
  lt_end_codons = ["TAA","TAG", "TGA"]
  lt = []
  for i in range(reading_frame, len(seq)-k+1, k):
    if seq[i:i+k] in lt_end_codons:
      lt.append((i, seq[i:i+k]))
  if lt:
    return lt
  else:
    return False


def read_seq_from_file(filename):
    """ Reads a sequence from a multi-line text file. """
    fh = open(filename, "r")
    lines = fh.readlines()
    seq = ""
    for l in lines:
        seq += l.replace("\n","")
    fh.close()
    return seq

def write_seq_to_file(seq, filename):
    """ Writes a sequence to file. """
    with open(filename, "w") as f:
      f.write(seq)
      print(f"{filename} written correctly")

def read_genetic_code_from_file(filename):
    """ Reads the genetic code to a dictionary from a multi-line text file. """
    import re
    with open(filename, "r") as f:
        file_content = f.read()
    seq = "".join(re.split("[^a-zA-Z]*", file_content))
    return seq

def validate_dna(dna_seq):
    """ Checks if DNA sequence is valid. Returns True is sequence is valid, or False otherwise. """
    seqm = dna_seq.upper()
    set_basis = {"A", "T", "C", "G"}
    s = (set(dna_seq) - set_basis)
    if (set_basis - set(dna_seq)):
        raise Exception(f"There are different characters in this sequence: {s}")
    else:
        return True

# def frequency(seq):
#     """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
#     dic = {}
#     for s in seq.upper():
#         if s in dic: dic[s] += 1
#         else: dic[s] = 1
#     return dic

def frequency(s, k, steps):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    d = {}
    for i in range(0, len(s) - k + 1, steps):
        if s[i:i+k] in list(d.keys()):
            d[s[i:i+k]]+=1
        else:
            d[s[i:i+k]]=1
    return d

def gc_content(dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def gc_content_subseq(dna_seq, k=100):
    """ Returns GC content of non-overlapping sub-sequences of size k. """
    # complete
    # ...


def transcription(dna_seq):
    """ Function that computes the RNA corresponding to the transcription of the DNA sequence provided. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    return dna_seq.upper().replace("T","U")


def reverse_complement(dna_seq, type_basis="dna"):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    if type_basis == "dna":
        dict_comp_basis = {"A":"T", "C":"G", "T":"A", "G":"C"}
    elif type_basis == "rna":
        dict_comp_basis = {"A":"U", "C":"G", "U":"A", "G":"C"}
    lt_comp = list(map(lambda x: dict_comp_basis[x],  dna_seq))
    return "".join(lt_comp[::-1])


def get_translate_codon_dict():
    return  {
      "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"
    }

def translate_codon(cod):
    """Translates a codon into an aminoacid using an internal dictionary with the standard genetic code."""
    tc = {"GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
      "TGT":"C", "TGC":"C",
      "GAT":"D", "GAC":"D",
      "GAA":"E", "GAG":"E",
      "TTT":"F", "TTC":"F",
      "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
      "CAT":"H", "CAC":"H",
      "ATA":"I", "ATT":"I", "ATC":"I",
      "AAA":"K", "AAG":"K",
      "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
      "ATG":"M", "AAT":"N", "AAC":"N",
      "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
      "CAA":"Q", "CAG":"Q",
      "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
      "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
      "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
      "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
      "TGG":"W",
      "TAT":"Y", "TAC":"Y",
      "TAA":"_", "TAG":"_", "TGA":"_"}
    if cod in tc: return tc[cod]
    else: return None


def translate_seq(dna_seq, reading_frame=0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    assert validate_dna(dna_seq)
    k = 3
    dict_codons_to_amino = get_translate_codon_dict()
    seq = dna_seq[reading_frame:]
    lt_seq_broken = [seq[i:i+k].upper() for i in range(0, len(seq) - k + 1, k)]
    map_codons_into_amino = map(lambda x: dict_codons_to_amino[x.replace("U","T")], lt_seq_broken)
    return "".join(list(map_codons_into_amino))


def codon_usage(dna_seq, aa):
    """Provides the frequency of each codon encoding a given aminoacid, in a DNA sequence ."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    seqm = dna_seq.upper()
    dic = {}
    total = 0
    for i in range(0, len(seqm)-2, 3):
        cod = seqm[i:i+3]
        if translate_codon(cod) == aa:
            if cod in dic:
                dic[cod] += 1
            else: dic[cod] = 1
            total += 1
    if total >0:
        for k in dic:
            dic[k] /= total
    return dic


def reading_frames(dna_seq):
  """Computes the six reading frames of a DNA sequence (including the reverse complement."""
  assert validate_dna(dna_seq)
  rc = reverse_complement(dna_seq, "dna")
  lt = [(translate_seq(dna_seq,i), translate_seq(rc,i)) for i in range(0,3)]

  return list(list(zip(*lt))[0] + list(zip(*lt))[1])


# def all_proteins_rf(aa_seq):
#     """Computes all posible proteins in an aminoacid sequence."""
#     aa_seq = aa_seq.upper()
#     current_prot = []
#     proteins = []
#     for aa in aa_seq:
#         if aa == "_":
#             if current_prot:
#                 for p in current_prot:
#                     proteins.append(p)
#                 current_prot = []
#         else:
#             if aa == "M":
#                 current_prot.append("")
#             for i in range(len(current_prot)):
#                 current_prot[i] += aa
#     return proteins

def all_proteins_rf(aa_seq, return_positions=False, max_len=False):
    """Computes all posible proteins in an aminoacid sequence."""
    aa_seq = aa_seq.upper()
    current_prot = []
    proteins = []
    positions = []
    current_pos = []
    for aa in range(0, len(aa_seq)):
        if aa_seq[aa] == "_":
            if current_prot:
                current_pos.append(aa)
                if len(current_prot) > 1 and max_len:
                    current_prot_temp = sorted(current_prot, key=lambda x: len(x[0]), reverse=True)
                    proteins.append(current_prot_temp[0])
                    positions.append([current_pos[0], current_pos[-1]])
                else:
                    for i in range(0, len(current_prot)):
                        proteins.append(current_prot[i])
                        positions.append([current_pos[i], current_pos[-1]])

                current_prot = []
                current_pos = []
        else:
            if aa_seq[aa] == "M":
                current_prot.append("")
                current_pos.append(aa)
            for i in range(len(current_prot)):
                current_prot[i] += aa_seq[aa]
                # current_pos.append(aa)
    # proteins = list(filter(lambda x: len(x)>=2, proteins))
    # positions = list(filter(lambda x: len(x)>=2, positions))
    if return_positions:
        return proteins, positions
    else: 
        return proteins

def all_orfs(dna_seq, return_positions=False):
    """Computes all possible proteins for all open reading frames."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    lt_all_rf = reading_frames(dna_seq)
    f = lambda x: all_proteins_rf(x, return_positions)
    if return_positions:
        lt_proteins = list(zip(*list(map(f, lt_all_rf))))[0]
        lt_positions = list(zip(*list(map(f, lt_all_rf))))[1]
        return sum(lt_proteins,[]), sum(lt_positions, [])
    else:
        return sum(list(map(f, lt_all_rf)),[])

def all_orfs_ord(dna_seq, minsize = 0):
    """Computes all possible proteins for all open reading frames. Returns ordered list of proteins with minimum size."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    res = all_orfs(dna_seq)
    res = list(filter(lambda x: len(x) >= minsize, res))
    res.sort(key=lambda s: len(s), reverse=True)
    return res

def insert_prot_ord(prot, list_prots):
    ''' inserts prot in list_prots in a sorted way '''
    i = 0
    while i < len(list_prots) and len(prot) < len(list_prots[i]):
        i += 1
    list_prots.insert(i, prot)


def test_frequency():
    seq_aa = input("Protein sequence:")
    freq_aa = frequency(seq_aa)
    list_f = sorted(freq_aa.items(), key=lambda x: x[1], reverse = True)
    for (k,v) in list_f:
        print("Aminoacid:", k, ":", v)

def test_all():
    seq = input("Insert DNA sequence:")
    if validate_dna (seq):
        print ("Valid sequence")
        print ("Transcription: ", transcription (seq))
        print("Reverse complement:", reverse_complement(seq))
        print("GC content (global):", gc_content(seq))
        print("Direct translation:" , translate_seq(seq))
        print("All proteins in ORFs (decreasing size): ", all_orfs_ord(seq))
    else: print("DNA sequence is not valid")

def test_files():
    fname = input("Insert input filename:")
    seq = read_seq_from_file(fname)
    if validate_dna (seq):
        print ("Valid sequence")
        print ("Transcription: ", transcription (seq))
        print("Reverse complement:", reverse_complement(seq))
        print("GC content (global):", gc_content(seq))
        print("Direct translation:" , translate_seq(seq))
        orfs = all_orfs_ord(seq)
        i = 1
        for orf in orfs:
            write_seq_to_file(orf, "orf-"+str(i)+".txt")
            i += 1
    else: print("DNA sequence is not valid")


if __name__ == "__main__":

    # test here all implemented functions
    # used your own defined sequences or read from example files
    seq = "ATGAGCGAC"
    print("Reverse:", seq[::-1])
    print("Reverse complement:", reverse_complement(seq))

    ## uncomment the test function to run
    test_frequency()
    test_all()
    #test_files()