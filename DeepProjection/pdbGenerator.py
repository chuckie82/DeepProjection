# SAXS simulation:
import os
import numpy as np
import matplotlib.pyplot as plt
import skopi as sk
import pypdb 

class PDBGenerator:
    def __init__(self):
       self.pdbID = None
       self.pdbFile = None

    def generate_pdb_id(self):
        # Randomly generate PDB ID (4 characters)
        # 1-9, [a-z,0-9], [a-z,0-9], [a-z,0-9]

        alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                 'n','o','p','q','r','s','t','u','v','w','x','y','z']
        num = ['0','1','2','3','4','5','6','7','8','9']
        alphanum = alpha + num
        self.pdbID = np.random.choice(num[1:]) + \
                     np.random.choice(alphanum) + \
                     np.random.choice(alphanum) + \
                     np.random.choice(alphanum)
        return self.pdbID

    def get_random_pdb(self):
        # Returns a random but valid PDB ID and file
        while 1:
            try:
                self.pdbID = self.generate_pdb_id()
                print("get_random_pdb trying:", self.pdbID)
                self.pdbFile = pypdb.get_pdb_file(self.pdbID, filetype='pdb', compression=False)
                return self.pdbID, self.pdbFile
            except:
                pass

    
    def download_pdb_in_range(self, lower, upper, dir):
        """
        Downloads all PDB files from the Protein Data Bank within the
        range of [lower, upper]

        Parameters
        ----------
        lower : str
            A PDB ID that will be the lower bound of our search.
            The lower range input will be considered a "wildcard" if the string contains
            an '*' after a value (e.g. 1A*) or the PDB ID is incomplete (e.g. 1A)
        
        upper : str
            A PDB ID that will be the upper bound of our search.
            The upper range input will be considered a "wildcard" if the string contains
            an '*' after a value (e.g. 1A*) or the PDB ID is incomplete (e.g. 1A)

        dir : str
            The directory to store the PDB files in.

        Corner Cases
        ------------

        Input Errors
        1. lower or upper are not of type str
        2. lower > upper; that is, lower comes after upper (e.g. BBBB to AAAA).
        3. The first character in the string for either lower or upper is a letter; not proper PDB ID format.
        4. len(lower) or len(upper) is greater than 4; not proper PDB ID format
        5. lower and upper contain upper case letters; just convert lower and upper to lower case.
        6. Keep in Mind: wildcard symbol or a space is found in between two sets of letters (e.g. 1A*B, *BBB, " AA" or "AD C")
        
        Special Ranges
        1. lower == upper; in this case, the function will attempt to download 1 file.
        2. lower == '*' and upper == '*'; in this case, we will download all files.
        3. lower is a valid ID, but upper == '*'; in this case, we will download all files between lower and upper == '9ZZZ'
        4. lower == '*', but upper is a valid ID; in this case, we will download all files between lower == '0000' and upper
        
        """
        
        alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                 'n','o','p','q','r','s','t','u','v','w','x','y','z']
        num = ['0','1','2','3','4','5','6','7','8','9']
        
        # Idea behind having the numbers go before the letters
        # is to count up from 0 to Z. So 0000 would be the first
        # ID, 0001 would be the second, and so on
        numalpha = num + alpha

        
        
        # Checking for Input Errors
        if isinstance(dir, str) == False:
            return False
        elif (isinstance(lower, str) == False) or (isinstance(upper, str) == False):    # Input Error 1: lower and upper are not str.
            return False
        else:
            # Input Error 5: lower and upper contain upper case letters; just convert lower and upper to lower case.
            lower = lower.lower()
            upper = upper.lower()


        if (lower > upper) and (lower != '*') and (upper != '*'):   # Input Error 2: lower > upper (MAY BE SOURCE OF BUGS)
            return False
        elif (lower[0] in alpha) or (upper[0] in alpha):            # Input Error 3: The first character in either lower or upper is a letter
            return False
        elif (len(lower) > 4) or (len(upper) > 4):                  # Input Error 4: len(lower) or len(upper) is greater than 4; not proper PDB ID format
            return False
        else:
            pass
        

        
        # Checking for Special Ranges
        # If lower == upper, our loop code should be able to handle that case on its own.

        if lower == '*' and upper == '*': # Special Range 2: lower == '*' and upper == '*', so we will download all files
            lower = '0000'
            upper = '9zzz'
        elif upper == '*':                # Special Range 3: lower is a valid ID, but upper == '*'
            upper = '9zzz'
        elif lower == '*':                # Special Range 4: lower == '*', but upper is a valid ID
            lower = '0000'
        else:
            pass

        
        
        # Now that the inputs were checked for Input Errors and Special Ranges,
        # we can preprocess the inputs for wildcards
        if '*' in lower:
            # If lower contains a '*', that means a wildcard is present.
            # Find the '*', and fill in the remaining characters with '0'
            for idx in range(len(lower)):
                if lower[idx] == '*':
                    lower = lower[0: idx] + (numalpha[0] * (4 - idx))
                    break
                else:
                    pass

        if len(lower) < 4:
            # If len(lower) < 4, ID contains a wildcard through incompletion.
            # Because length is less than 4, we will add the remaining characters as '0'.
            lower = lower + (numalpha[0] * (4 - len(lower)))

        if '*' in upper:
            # If upper contains a '*', that means a wildcard is present.
            # Find the '*', and fill in the remaining characters with '0'
            for idx in range(len(upper)):
                if upper[idx] == '*':
                    upper = upper[0: idx] + (numalpha[0] * (4 - idx))
                    break
                else:
                    pass
        
        if len(upper) < 4:
            # If len(upper) < 4, ID contains a wildcard through incompletion.
            # Because length is less than 4, we will add the remaining characters as '0'.
            upper = upper + (numalpha[0] * (4 - len(upper)))
        
        
        # Inputs have be preprocessed. Time to get the PDB files
        lower_firstpos = numalpha.index(lower[0])
        lower_secondpos = numalpha.index(lower[1])
        lower_thirdpos = numalpha.index(lower[2])
        lower_fourthpos = numalpha.index(lower[3])

        currentID = ''

        for firstpos in range (lower_firstpos, 10): # Have to add the +1 to the range input, as we want upper to be included in search
            
            if currentID > upper:
                # currentID is greater than upper, stop loop
                break   
            
            for secondpos in range(lower_secondpos, 36):
                
                if currentID > upper:
                    # currentID is greater than upper, stop loop
                    break   
                
                for thirdpos in range(lower_thirdpos, 36):
                    
                    if currentID > upper:
                        # currentID is greater than upper, stop loop
                        break   
                    
                    for fourthpos in range(lower_fourthpos, 36):
                        
                        currentID = numalpha[firstpos] + numalpha[secondpos] + numalpha[thirdpos] + numalpha[fourthpos]

                        if currentID > upper:
                            # currentID is greater than upper, stop loop
                            break  
                        else:

                            # Try to find the PDB file with ID of currentID and write to disk
                            try:
                                print("download_pdb_in_range trying:", currentID)
                                
                                pdb = pypdb.get_pdb_file(currentID, filetype='pdb', compression=False)
                                
                                with open(os.path.join(dir, currentID + ".pdb") ,"w") as f:
                                    f.writelines(pdb)

                                print("downloaded: ", currentID)
                            except:
                                pass

        # Return True to state that operation was a success
        return True



