import os
import numpy as np
import matplotlib.pyplot as plt
import pdbGenerator as pg
import pypdb
import skopi as sk

showPlot = 0
maxPDBs = 30000
outdir = "/reg/data/ana03/scratch/yoon82/Software/DeepProjection/DeepProjection/data"

Pdb = pg.PDBGenerator()

numPDBs = 0
while numPDBs < maxPDBs:
    # Search PDB by structure similarity
    found_pdbs = None
    while not found_pdbs:
        randPDB, pfile = Pdb.get_random_pdb()
        if not os.path.exists(os.path.join(outdir,randPDB+".npz")):
            try:
                found_pdbs = pypdb.Query(randPDB, query_type="structure").search()
            except:
                pass
    numNeighbors = 8
    if len(found_pdbs) < numNeighbors:
        numNeighbors = len(found_pdbs)

    # Fetch pdb file
    for i,val in enumerate(found_pdbs[:numNeighbors]):
        if not os.path.exists(os.path.join(outdir,val+".npz")):
            print("fetching: ", val)
            try:
                pdb_file = pypdb.get_pdb_file(val, filetype='pdb', compression=False)
                with open(os.path.join(outdir,val+".pdb"),"w") as f:
                    f.writelines(pdb_file)
            except:
                print("Could not fetch pdb file: ", val)
                continue

        # Set up particle
        particle = sk.Particle()
        if os.path.exists(os.path.join(outdir,val+".pdb")):
            try:
                particle.read_pdb(os.path.join(outdir,val+".pdb"), ff='WK')
            except:
                continue
        else:
            continue
        N = 100000 # no. of random HKL samples
        resmax = 1e-9 # maximum resolution of SAXS curve (m)
        saxs = sk.SAXS(particle,N,resmax)
        np.savez(os.path.join(outdir,val+".npz"),qs=saxs.qs,saxs=saxs.saxs,qmax=saxs.qmax)
        os.remove(os.path.join(outdir,val+".pdb"))
        numPDBs += 1

        if showPlot:
            plt.subplot(2,4,i+1)
            plt.yscale('log')
            plt.xlim(0,saxs.qmax/10**10) # convert to Angstroem
            plt.xlabel('q (inverse Angstroem)')
            plt.ylabel('logI')
            plt.plot(saxs.qs/10**10, saxs.saxs) # convert to Angstroem
            plt.title(val)

    if showPlot:
        plt.subplots_adjust(hspace=0.4, wspace=0.4) # spread out plots
        plt.show()


