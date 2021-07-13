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

def test_download_pdb_in_range(empty_pdbgenerator):
    '''download_pdb_in_range() can be run using an "empty" PDBGenerator'''
    '''download_pdb_in_range() returns True when it encounters no input errors; does not mean it successfully downloaded a file'''


    """
    Default values
        1. If upper == None, then range is from [lower, 9ZZZ].
        2. If dir == None, then download PDB files to current working directory.
    """
    # With lower == '0aaa', we download_pdb_in_range would equal False. That is because the first
    # character is a number from 1 to 9.
    # These tests look at Input Error #5 and Input Error #6
    assert empty_pdbgenerator.download_pdb_in_range('0aaa', None, None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0AAA', None, None) == False


    """
    Inputs with complete IDs; no wildcards
    """
    # Complete PDB IDs.
    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bbb', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bbb', 'FAKEDIR') == False

    # PDB IDs with wildcard character '*'.
    assert empty_pdbgenerator.download_pdb_in_range('0aa*', '0bbb', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aa*', '0bbb', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bb*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bb*', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('0aa*', '0bb*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aa*', '0bb*', 'FAKEDIR') == False

    # Incomplete PDB IDs, which are treated as "wildcards" by the function.
    assert empty_pdbgenerator.download_pdb_in_range('0aa', '0bbb', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aa', '0bbb', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bb', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0aaa', '0bb', 'FAKEDIR') == False


    """
    Input Errors
        1. lower or upper are not of type str
        2. lower > upper; that is, lower comes after upper (e.g. BBBB to AAAA).
        3. The first character in the string for either lower or upper is a letter; not proper PDB ID format.
        4. len(lower) or len(upper) is greater than 4; not proper PDB ID format.
        5. lower and upper contain upper case letters; just convert lower and upper to lower case.
        6. The first character is '0'; not proper PDB format.
    """

    # Testing Input Error #1: lower or upper are not of type str
    assert empty_pdbgenerator.download_pdb_in_range('1*', 2, None) == False
    assert empty_pdbgenerator.download_pdb_in_range(1, '2*', None) == False


    # Testing Input Error #2: lower > upper; that is, lower comes after upper (e.g. 2BBB to 1AAA)
    assert empty_pdbgenerator.download_pdb_in_range('2BBB', '1AAA', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('2BBB', '1AAA', 'FAKEDIR') == False

    # Testing Input Error #3: The first character in the string for either lower or upper is a letter; not proper PDB ID format.
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', 'BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', 'BBBB', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1AAA', 'BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAA', 'BBBB', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('AAAA', '1BBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', '1BBB', 'FAKEDIR') == False

    # Testing Input Error #4: len(lower) or len(upper) is greater than 4; not proper PDB ID format.
    assert empty_pdbgenerator.download_pdb_in_range('1AAAA', '2BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAAA', '2BBBB', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1AAA*', '2BBB*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAA*', '2BBB*', 'FAKEDIR') == False
