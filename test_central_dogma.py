from central_dogma import *


def test_antisense():
    assert antisense("TTAGGGCATG") == "CATGCCCTAA"


def test_transcribe_to_rna():
    assert transcribe_to_rna("TTAGGGCATG") == 'CAUGCCCUAA'


def test_codon_to_prot():
    assert codon_to_prot("AUG") == "Met"

def test_prot3_to_prot1():
    assert prot3_to_prot1("Met")=="M"


def test_translation():
    assert translation('CAUGCCCUAA') == ("MetPro",'MP')