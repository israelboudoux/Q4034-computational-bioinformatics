import pandas as pd
import re
import difflib
import sys

folder_path = "./"

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

def translate_seq(dna_seq, reading_frame=0):
    """ Translates a DNA sequence into an aminoacid sequence. """
    assert validate_dna(dna_seq)
    k = 3
    dict_codons_to_amino = get_translate_codon_dict()
    seq = dna_seq[reading_frame:]
    lt_seq_broken = [seq[i:i+k].upper() for i in range(0, len(seq) - k + 1, k)]
    map_codons_into_amino = map(lambda x: dict_codons_to_amino[x.replace("U","T")], lt_seq_broken)
    return "".join(list(map_codons_into_amino))

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

def reading_frames(dna_seq):
  """Computes the six reading frames of a DNA sequence (including the reverse complement."""
  assert validate_dna(dna_seq)
  rc = reverse_complement(dna_seq, "dna")
  lt = [(translate_seq(dna_seq,i), translate_seq(rc,i)) for i in range(0,3)]

  return list(list(zip(*lt))[0] + list(zip(*lt))[1])

def reverse_complement(dna_seq, type_basis="dna"):
    """ Computes the reverse complement of the inputted DNA sequence. """
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    if type_basis == "dna":
        dict_comp_basis = {"A":"T", "C":"G", "T":"A", "G":"C"}
    elif type_basis == "rna":
        dict_comp_basis = {"A":"U", "C":"G", "U":"A", "G":"C"}
    lt_comp = list(map(lambda x: dict_comp_basis[x],  dna_seq))
    return "".join(lt_comp[::-1])

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

def gc_content(dna_seq):
    """ Returns the percentage of G and C nucleotides in a DNA sequence. """
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc": gc_count += 1
    return gc_count / len(dna_seq)

def frequency(s, k, steps):
    """ Calculates the frequency of each symbol in the sequence. Returns a dictionary. """
    d = {}
    for i in range(0, len(s) - k + 1, steps):
        if s[i:i+k] in list(d.keys()):
            d[s[i:i+k]]+=1
        else:
            d[s[i:i+k]]=1
    return d

def get_orf_names(lt_coords):
    import string
    lt_letters = string.ascii_letters[:26]
   
    current_position = []
    ORFs = []
    i = 1
    j = 1
    k = 1

    count = 0

    for start, end in lt_coords:
        count+=1
        if current_position:
            # Case 1: ORF inside ORF
            if start >= current_position[0] and end <= current_position[1]:
                # ORFs.append(f"ORF{i}.{j}, ({start}, {end})")
                ORFs.append(f"ORF{i}.{j}")
                j+=1
            # Case 2: ORF overlaps ORF
            elif start < current_position[1] and end > current_position[1]:
                letter = lt_letters[k]
                # ORFs.append(f"ORF{i}.{letter}, ({start}, {end})")
                ORFs.append(f"ORF{i}.{letter}")
                k+=1
            # Case 3: ORF different ORF
            elif start > current_position[1]:
                j = 1
                k = 1
                i+=1

                current_position = []
                current_position.append(start)
                current_position.append(end)

                # ORFs.append(f"ORF{i}, ({start}, {end})")
                ORFs.append(f"ORF{i}")
        else:
            current_position.append(start)
            current_position.append(end)

            # ORFs.append(f"ORF{i} ({start}, {end})")
            ORFs.append(f"ORF{i}")

    if len(ORFs)!=count:
        raise Exception("Some ORFs weren't counted.")
    else:
        return ORFs

def validate_dna(dna_seq):
  seqm = dna_seq.upper()
  set_basis = {"A", "T", "C", "G"}
  s = (set(dna_seq) - set_basis)
  if (set_basis - set(dna_seq)):
      raise Exception(f"There are different characters in this sequence: {s}")
  else:
      return True

def read_fasta_2dictionary(filename):
  with open(filename, "r") as f:
    f = f.read()

  import re
  lt_idx = re.findall(r"\>(.*?)\ ",f)

  lt_entries = f.split(">")
  lt_seqs = []
  for i in range(0, len(lt_entries)):
    if lt_entries[i]:
      raw_seq = lt_entries[i][lt_entries[i].find("\n"):]
      clean_seq = "".join(re.split("[^a-zA-Z]*", raw_seq)).upper()
      assert validate_dna(clean_seq)
      lt_seqs.append(clean_seq)

  return dict(zip(lt_idx, lt_seqs))

def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    return s1[pos_a:pos_a+size]

def main(sequence):
    pdf_proteins = pd.read_csv(folder_path+"proteins_86693_757732.csv")

    '''
    =====================================//=====================================
    STEP 00 - Get statistics
    =====================================//=====================================
    In this step, we are applying all the filters stipulated by the proposed 
    exercise, namely:
    ------------------------------------//--------------------------------------
    '''
    # 1) Length of the sequence.
    print(f"\n -> Sequence len: {len(sequence)}")

    # 2) Frequency in %:
    print(f"\n -> Frequence:")
    freq = frequency(sequence, 1, 1)
    total = sum(list(freq.values()))
    for i, j in freq.items():
        print(f"    {i}: {str(round(j/total*100, 2))}%")

    # 3) GC content:
    print(f"\n -> GC content: {round(gc_content(sequence)*100, 2)}%")

    # 4,5) Number of start/end codon:
    print("\n -> Number of start/end codons:")
    for i in range(0,3):
        l_start = len(start_codon_position(sequence, i))
        l_end = len(start_codon_position(reverse_complement(sequence), i))
        print(f"    * In frame {i}: \n      - Start: {l_start} | Stop: {l_end}")

    # 6) Most and less freq. codon:
    dict_freq_codons = frequency(sequence, 3, 1)
    lt = list(dict_freq_codons.values())
    lt_max = max(lt)
    lt_min = min(lt)

    print("\n -> Most and less freq. codon:")
    for i, j in  dict_freq_codons.items():
        if j == lt_max or j == lt_min:
            print(f"    {i}: {j}")
    
    ''' 
    =====================================//=====================================
    STEP 01 - Creating main dataset
    =====================================//=====================================
    In this step, we are creating a dataframe that consolidates all the information 
    regarding the ORFs and their respective positions, for each of the possible DNA 
    sequences (3 direct frames and 3 reverse frames).
    ------------------------------------//--------------------------------------
    '''
    #  Getting 6 sequences, one for each direct and reverse reading frame:
    lt_reading_frames = reading_frames(sequence)
    # Putting those sequences into a Pandas DataFrame with "reading_frame" as index:
    pdf = pd.DataFrame(lt_reading_frames, columns=["seq__amino"])
    pdf.index.names = ["reading_frame"]
    # Getting all ORFs associated with those 6 sequences with "all_proteins_rf" function;
    # Creating two new columns with the outputs of "all_proteins_rf" (orfs and positions);
    # Getting only the ORFs that have the longest length when ORF overlap occurs (max_len):
    pdf["seq__orf"], pdf["orf__coord"] = zip(*pdf.apply(lambda x: all_proteins_rf(x["seq__amino"], True, max_len=False), axis=1))
    # Exploding the proteins and positions calculated above into rows:
    pdf_exp = pdf.explode(["seq__orf", "orf__coord"])
    pdf_exp["orf__coord_start"], pdf_exp["orf__coord_end"] = zip(*pdf_exp["orf__coord"])
    # Creating a column with the length of each ORF:
    pdf_exp["seq__orf_len"] = pdf_exp["seq__orf"].str.len()
    # Creating two new columns for dna coordinates calculated from orf/aminoacid coordinates:
    for n in [0,1,2]:
        pdf_exp.loc[n, "dna__coord_start"] = pdf_exp.loc[n, "orf__coord_start"]*3 + n
        pdf_exp.loc[n, "dna__coord_end"] = pdf_exp.loc[n, "orf__coord_end"]*3 + 3 + n
    #droping reading frames of the reverse sequence
    pdf_exp = pdf_exp.dropna()
    # Getting the DNA sub-sequences associated with start and end of each ORF:
    pdf_exp["seq__dna"] = pdf_exp.apply(lambda x: sequence[int(x["dna__coord_start"]):int(x["dna__coord_end"])], axis=1)

    '''
    =====================================//=====================================
    STEP 02 - Applying filters
    =====================================//=====================================
    In this step, we are applying all the filters stipulated by the proposed 
    exercise, namely:
    (1) sequences must contain >= 40 amino acids;
    (2) If an ORF is contained within another, select only the largest of them;
    (3) order the data by ORF length.
    ------------------------------------//--------------------------------------
    '''
    # Creating a new dataframe "prf_exp_temp" to apply filters and other particular manipulations:
    pdf_exp_temp = pdf_exp.copy()
    # Filtering only ORFs with at least 40 aminoacis in its sequence:
    pdf_exp_temp = pdf_exp_temp[pdf_exp_temp["seq__orf"].str.len() >= 40]
    # Ordering dataframe by ORF coordinates to apply the function that will create "names"/"labes" for those ORFs:
    names_by_frame = True
    if names_by_frame:
        pdf_exp_temp = pdf_exp_temp.reset_index().sort_values(["reading_frame","orf__coord_start","orf__coord_end"], ascending=[True,True,False]).set_index("reading_frame")
        lt_orf_names = pdf_exp_temp.groupby("reading_frame").apply(lambda x: get_orf_names(x[["dna__coord_start", "dna__coord_end"]].values))
        lt_orf_names = [[f"rf0{j}_" + i for i in lt_orf_names[j]] for j in range(0, len(lt_orf_names))]
        pdf_exp_temp["orf__name"] = sum(lt_orf_names, [])
    else:
        pdf_exp_temp = pdf_exp_temp.sort_values("orf__coord")
        pdf_exp_temp["orf__name"] = get_orf_names(pdf_exp_temp[["dna__coord_start", "dna__coord_end"]].values)
    # Sorting dataframe by orf length:
    pdf_exp_temp = pdf_exp_temp.sort_values("seq__orf_len", ascending=False)
    # Just sorting columns alphabeticaly:
    pdf_exp_temp = pdf_exp_temp[sorted(pdf_exp_temp.columns)]
    
    print("\n -> creating all_potential_proteins.txt file...")
    (
        pdf_exp_temp
            .sort_values("seq__orf_len", ascending=False)["seq__dna"]
            .to_csv(folder_path+"all_potential_proteins.txt", header=False, index=False)
    )
    print("    * file created!")

    print("\n -> creating orf_coordinates.txt file...")
    
    pattern = "rf0\d{1}_ORF\d+"
    l = pd.Series([re.findall(pattern, orf)[0] for orf in pdf_exp_temp["orf__name"] if re.findall(pattern, orf)]).unique()
    pdf_exp_temp_ex08 = pdf_exp_temp[pdf_exp_temp["orf__name"].isin(l)].sort_values("orf__coord")
    pdf_exp_temp_ex08 = pdf_exp_temp_ex08.sort_values("seq__orf_len", ascending=False)
    (
        pdf_exp_temp_ex08[["dna__coord_start", "dna__coord_end", "seq__dna"]]
            .to_csv(
                folder_path+"orf_coordinates.txt", 
                float_format="%.0f", 
                sep=",", 
                index=False, 
                header=False
            )
    )
    print("    * file created!")

    '''
    =====================================//=====================================
    STEP 03 - Applying filters
    =====================================//=====================================
    In this step, we are applying all the filters stipulated by the proposed 
    exercise, namely:
    (1) sequences must contain >= 40 amino acids;
    (2) If an ORF is contained within another, select only the largest of them;
    (3) order the data by ORF length.
    ------------------------------------//--------------------------------------
    '''
    pdf_proteins["seq__dna"] = pdf_proteins.apply(lambda x: sequence[int(x["Start"])-1:int(x["Stop"])], axis=1)
    pdf_proteins["dna__coord_start"] = pdf_proteins["Start"] - 1
    pdf_proteins["dna__coord_end"] = pdf_proteins["Stop"]
    pdf_temp = (
        pdf_proteins[
            ["Locus", "dna__coord_start", "dna__coord_end", "seq__dna"]
            ]
        .merge(
            pdf_exp_temp[
                ["seq__dna", "dna__coord_start", "dna__coord_end", "orf__name"]
                ], 
            how="cross"
            )
        )

    print("\n -> Getting overlap between ORFs...")
    pdf_temp["overlap__seq"] = pdf_temp.apply(lambda x: get_overlap(str(x["seq__dna_x"]), str(x["seq__dna_y"])), axis=1)
    pdf_temp["overlap__perc"] = pdf_temp["overlap__seq"].str.len()/pdf_temp["seq__dna_x"].str.len()
    print("    * Overlaps calculated!")
    for i in pdf_temp.sort_values(["Locus", "overlap__perc"], ascending=[True, False]).groupby("Locus").head(1)[["Locus","overlap__perc"]].values:
        print(f"      - Locus: {i[0]} > Overlap: {i[1]}")

if __name__ == "__main__":
    try:
        filename = str(sys.argv[1])
    except:
        raise Exception("You need to pass <filename.extension> as a parameter!")

    dict_sequence = read_fasta_2dictionary(folder_path+filename)
    lt_keys = list(enumerate(list(dict_sequence.keys())))
    msg = f"""
    Select one of the numbers below by typing and then pressing Enter:
    {lt_keys}
    """
    num = int(input(msg))
    key = list(filter(lambda x: x[0] == num, lt_keys))[0][1]
    sequence = dict_sequence[key]
    main(sequence)