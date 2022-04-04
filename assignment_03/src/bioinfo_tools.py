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

def number_of_proteins_by_frame(seq):
  lt_rf = bsf.reading_frames(seq)
  lt_prots = list(map(bsf.all_proteins_rf, lt_rf))

  print("-> Number of possible proteins:")
  print("  > direct sequence")
  j = 0
  for i in range(0, len(lt_prots)):
    print(f"    * In frame {j}: {len(lt_prots[i])}")
    j+=1
    if i == 2: 
      j=0
      print("  > reverse sequence")

def get_orf_names(lt_coords):
    import string
    lt_letters = string.ascii_letters[:26]
    """
    CASE 3:
    ORF1 ======== *********** ORF2

    CASE 2:
    ORF1a =========
               ********* ORF1b
    
    CASE 1:
    ORF1ab ================
             ********* ORF1b
    """
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