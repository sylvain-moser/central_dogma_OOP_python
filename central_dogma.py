
def antisense(seq):
  translate = {"A":"T", "C":"G", "G":"C", "T":"A"}
  rev = seq[::-1]
  return ''.join([translate[letter] for letter in rev ])
  # [k, v for k, v in translate]

def translate_to_rna(seq):
  translate = {"A":"U", "C":"G", "G":"C", "T":"A"}
  rev = seq[::-1]
  return ''.join([translate[letter] for letter in rev ])

def all_rna(seq):
  return [translate_to_rna(seq), translate_to_rna(antisense(seq))]

class RNAseq:
    def __init__(self,sequence,sense):
        self.sequence=sequence
        self.sense=sense
    def transcribe(self):
        translate = {"A": "U", "C": "G", "G": "C", "T": "A"}
        if self.sense=="+":
            rev = self.sequence[::-1]
            return ''.join([translate[letter] for letter in rev])
        elif self.sense=="-":
            sense_to_antisens = {"A": "T", "C": "G", "G": "C", "T": "A"}
            rev = self.sequence[::-1]
            other_sense=''.join([sense_to_antisens[letter] for letter in rev])
            return ''.join([translate[letter] for letter in other_sense])




seq=RNAseq("ATCGCTA","+")
trans=seq.transcribe()
print(trans)
seq2=seq=RNAseq("ATCGCTA","-")
trans2=seq2.transcribe()
print (trans2)

