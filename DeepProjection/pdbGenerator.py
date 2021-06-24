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



