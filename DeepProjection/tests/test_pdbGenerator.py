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
    # There is no PDB with an ID of '9zzz', so it will be used in default value test.
    assert empty_pdbgenerator.download_pdb_in_range('9zzz', None, None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZZ', None, None) == True

    """
    Starting off with some tests on inputs with no errors
    '9zzy' and '9zzz' do not exist in Protein Data Bank, perfect for testing.
    """
    # Complete PDB IDs.
    assert empty_pdbgenerator.download_pdb_in_range('9zz0', '9zzz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zz0', '9ZZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zz0', '9Zzz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zz0', '9zZZ', None) == True

    assert empty_pdbgenerator.download_pdb_in_range('9zzy', '9zzz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZY', '9ZZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zzy', '9zzz', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZY', '9ZZZ', 'FAKEDIR') == True
    
    # PDB IDs with wildcard character '*'.
    assert empty_pdbgenerator.download_pdb_in_range('9zz*', '9zzz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZ*', '9ZZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zz*', '9zzz', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZ*', '9ZZZ', 'FAKEDIR') == True

    assert empty_pdbgenerator.download_pdb_in_range('9zyz', '9zz*', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zyz', '9ZZ*', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zyz', '9zz*', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zyz', '9ZZ*', 'FAKEDIR') == True

    assert empty_pdbgenerator.download_pdb_in_range('9z*', '9z1*', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z*', '9Z1*', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9z*', '9z1*', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z*', '9Z1*', 'FAKEDIR') == True

    # Incomplete PDB IDs, which are treated as "wildcards" by the function.
    assert empty_pdbgenerator.download_pdb_in_range('9zz', '9zzz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZ', '9ZZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zz', '9zzz', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9ZZ', '9ZZZ', 'FAKEDIR') == True

    assert empty_pdbgenerator.download_pdb_in_range('9zyz', '9zz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zyz', '9ZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9zyz', '9zz', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zyz', '9ZZ', 'FAKEDIR') == True

    # These lines are failing because '9z00' is > '9z', so it triggers one of the error checks.
    # Best thing to do is have the wildcards be processed first, and then checked afterwards.
    assert empty_pdbgenerator.download_pdb_in_range('9z00', '9z', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z00', '9Z', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9z00', '9z', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z00', '9Z', 'FAKEDIR') == True
    ################################################################################


    assert empty_pdbgenerator.download_pdb_in_range('9z', '9zz', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z', '9ZZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9z', '9zz', 'FAKEDIR') == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z', '9ZZ', 'FAKEDIR') == True

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
    assert empty_pdbgenerator.download_pdb_in_range('1*', 2, 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range(1, '2*', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1', 2, None) == False
    assert empty_pdbgenerator.download_pdb_in_range(1, '2', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1', 2, 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range(1, '2', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1000', 2000, None) == False
    assert empty_pdbgenerator.download_pdb_in_range(1000, '2000', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1000', 2000, 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range(1000, '2000', 'FAKEDIR') == False


    # Testing Input Error #2: lower > upper; that is, lower comes after upper (e.g. 2BBB to 1AAA)
    assert empty_pdbgenerator.download_pdb_in_range('2BBB', '1AAA', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('2BBB', '1AAA', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('2BB*', '1AA*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('2BB*', '1AA*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('2BB', '1AA', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('2BB', '1AA', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('2000', '1000', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('2000', '1000', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('200*', '100*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('200*', '100*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('200', '100', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('200', '100', 'FAKEDIR') == False

    # Testing Input Error #3: The first character in the string for either lower or upper is a letter; not proper PDB ID format.
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', 'BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', 'BBBB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA*', 'BBB*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA*', 'BBB*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA', 'BBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA', 'BBB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('A*', 'B*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('A*', 'B*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('A', 'B', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('A', 'B', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1AAA', 'BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAA', 'BBBB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1AA*', 'BBB*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AA*', 'BBB*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1AA', 'BBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AA', 'BBB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1*', 'B*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1*', 'B*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1', 'B', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1', 'B', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('AAAA', '1BBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAAA', '1BBB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA*', '1BB*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA*', '1BB*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA', '1BB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('AAA', '1BB', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('A*', '1*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('A*', '1*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('A', '1', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('A', '1', 'FAKEDIR') == False

    # Testing Input Error #4: len(lower) or len(upper) is greater than 4; not proper PDB ID format.
    assert empty_pdbgenerator.download_pdb_in_range('10000', '20000', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('10000', '20000', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAAA', '2BBBB', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAAA', '2BBBB', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('1000*', '2000*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1000*', '2000*', 'FAKEDIR') == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAA*', '2BBB*', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('1AAA*', '2BBB*', 'FAKEDIR') == False

    # Testing Input Error #5: lower and upper contain upper case letters; just convert lower and upper to lower case.
    # These download_pdb_in_range() calls should return True, as the code should convert all strings to lowercase.
    assert empty_pdbgenerator.download_pdb_in_range('9Z00', '9Z01', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Z00', '9Z01', 'FAKEDIR') == True

    assert empty_pdbgenerator.download_pdb_in_range('9Zz0', '9Zz1', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9Zz0', '9Zz1', 'FAKEDIR') == True

    assert empty_pdbgenerator.download_pdb_in_range('9AAY', '9AAZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9AAY', '9AAZ', 'FAKEDIR') == True
    
    assert empty_pdbgenerator.download_pdb_in_range('9AaY', '9AaZ', None) == True
    assert empty_pdbgenerator.download_pdb_in_range('9AaY', '9AaZ', 'FAKEDIR') == True

    # Testing Input Error #6: The first character is '0'; not proper PDB format.
    assert empty_pdbgenerator.download_pdb_in_range('0Z00', '0Z01', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0Z00', '0Z01', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('7Zz0', '0Zz1', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('7Zz0', '0Zz1', 'FAKEDIR') == False

    assert empty_pdbgenerator.download_pdb_in_range('0AAY', '7AAZ', None) == False
    assert empty_pdbgenerator.download_pdb_in_range('0AAY', '7AAZ', 'FAKEDIR') == False
