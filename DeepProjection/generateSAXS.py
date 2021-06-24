import os
import matplotlib.pyplot as plt
import pdbGenerator as pg
import pypdb
import skopi as sk

Pdb = pg.PDBHandler()

# Search PDB by structure similarity
found_pdbs = None
while not found_pdbs:
    randPDB, pfile = Pdb.get_random_pdb()
    print("Got pdb: ", randPDB)
    try:
        print("Query")
        found_pdbs = pypdb.Query(randPDB, query_type="structure").search()
        print("found pdbs: ", found_pdbs)
    except:
        pass
numNeighbors = 8
if len(found_pdbs) < 8:
    numNeighbors = len(found_pdbs)

# Fetch pdb file
for i,val in enumerate(found_pdbs[:numNeighbors]):
    if not os.path.exists(val+".pdb"):
        print("fetching: ", val)
        pdb_file = pypdb.get_pdb_file(val, filetype='pdb', compression=False)
        with open(val+".pdb","w") as f:
            f.writelines(pdb_file)

    # Set up particle
    particle = sk.Particle()
    particle.read_pdb(val+".pdb", ff='WK')
    N = 100000 # no. of random HKL samples
    resmax = 1e-9 # maximum resolution of SAXS curve (m)
    saxs = sk.SAXS(particle,N,resmax)

    plt.subplot(2,4,i+1)
    plt.yscale('log')
    plt.xlim(0,saxs.qmax/10**10) # convert to Angstroem
    plt.xlabel('q (inverse Angstroem)')
    plt.ylabel('logI')
    plt.plot(saxs.qs/10**10, saxs.saxs) # convert to Angstroem
    plt.title(val)
plt.subplots_adjust(hspace=0.4, wspace=0.4) # spread out plots
plt.show()


