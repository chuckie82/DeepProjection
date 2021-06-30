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

    
    def download_all_existing_pdb(self):
        # Downloads all PDB files from the Protein Data Bank

        alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                 'n','o','p','q','r','s','t','u','v','w','x','y','z']
        num = ['0','1','2','3','4','5','6','7','8','9']
        
        # Idea behind having the number go before the letters
        # is to have 0000 be the starting ID, and go on from there.
        numalpha = num + alpha
        
        list_of_valid_ids = []
        currentID = '0000'

        # Time to generate the codes
        # The first position is always a number from 0 - 9
        for first_position in range (10):

            # The second position and all consecutive positions are a mix of letters and numbers
            for second_position in range(36):
                for third_position in range(36):
                    for forth_position in range(36):
                        currentID = numalpha[first_position] + numalpha[second_position] + numalpha[third_position] + numalpha[forth_position]
                        print(currentID)
                        list_of_valid_ids.append(currentID)

        # For debug
        print(len(list_of_valid_ids))

        return 0



