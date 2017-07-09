# README #

Python script that remove extra waters from an AMBER protein-membrane system. Tested with python 2.7. It assumes a PDB file with TER cards.

Examples are provided in the script remove_waters.py, such as:

./remove_waters.py -p proteinmemb.pdb -o out -s "resname LA PC" 

Extra features will be added soon.

Various python scripts to calculate PCA, Allosteric Networks, RMSD, and RMSF plots. It uses MDtraj, Networkx, numpy and scikit-learn.

Jose Flores-Canales