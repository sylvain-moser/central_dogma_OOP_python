import json

def antisense(seq):
  translate = {"A":"T", "C":"G", "G":"C", "T":"A"}
  rev = seq[::-1]
  return ''.join([translate[letter] for letter in rev ])
  # [k, v for k, v in translate]

def transcribe_to_rna(seq):
  translate = {"A":"U", "C":"G", "G":"C", "T":"A"}
  rev = seq[::-1]
  return ''.join([translate[letter] for letter in rev ])


def read_json_dictionary(jfile):
    file = open(jfile)
    dict = json.load(file)
    file.close()
    return (dict)


cod2prot = read_json_dictionary('codons.json')
prot2p1 = read_json_dictionary('peptides.json')


def codon_to_prot(codon):
    cod2prot = read_json_dictionary('codons.json')
    return cod2prot[codon]


def prot3_to_prot1(p3):
    prot2p1 = read_json_dictionary('peptides.json')
    return prot2p1[p3.lower()]


def translation(seq):
  start_seq = "AUG"             # prior knowledge
  start_pos = str.index(seq, start_seq)
  if start_pos == -1:
    return False
  peptide = codon_to_prot(start_seq)
  peptide1=prot3_to_prot1(peptide)
  for codon_pos in range(start_pos+3,    # codon after found start
                         len(seq)-2,     # last possible codon pos
                         3):             # in steps of 3
    protein = codon_to_prot(seq[codon_pos:codon_pos+3])
    if protein == 'Stop':          # STOP
      return peptide,peptide1
    protein1 = prot3_to_prot1(protein)
    peptide += protein
    peptide1 +=protein1
  return peptide,peptide1              # rather fail? (see Specification 3e)

class DNAseq:
    def __init__(self,sequence,sense):
        self.sequence=sequence
        self.sense=sense
    def transcribe(self):
        if self.sense=="+":
            return RNAseq(sequence=transcribe_to_rna(self.sequence))
        elif self.sense=="-":
            return RNAseq(sequence=transcribe_to_rna(antisense(self.sequence)))

class RNAseq:
    def __init__(self,sequence):
        self.sequence=sequence
    def translate(self):
        return Peptide(aa_sequence=translation(self.sequence)[0],short_aa_sequence=translation(self.sequence)[1])

class Peptide:
    def __init__(self,aa_sequence,short_aa_sequence):
        self.aa_sequence=aa_sequence
        self.short_aa_sequence=short_aa_sequence



seq=DNAseq("AGGACGGGCTAACTCCGCTCGTCACAAAGCGCAATGCAGCTATGGCAGATGTTCATGCCG","+")
trans=seq.transcribe()
print(trans.sequence)
prot1=trans.translate()
print(prot1.aa_sequence)
print(prot1.short_aa_sequence)

