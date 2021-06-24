import pytest
from pdbGenerator import PDBGenerator

@pytest.fixture
def empty_pdbgenerator():
    '''Returns a PDBGenerator instance with no pdb'''
    return PDBGenerator()

def test_default_initial_pdbgenerator(empty_pdbgenerator):
    assert empty_pdbgenerator.pdbID == None

def test_getting_pdb(empty_pdbgenerator):
    pdbID, pdbFile = empty_pdbgenerator.get_random_pdb()
    assert empty_pdbgenerator.pdbID is not None


