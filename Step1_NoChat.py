from Bio.PDB import PDBParser
from Bio.PDB import *

def get_interface(structure,threshold):
    chainA=structure[0]["A"]
    chainE=structure[0]["E"]
    residueSetA=set()
    residueSetE=set()
    for i in chainA.get_residues():
        for y in i.get_atoms():
            for j in chainE.get_residues():
                for z in j.get_atoms():                
                    if y-z <= threshold:
                        residueSetA.add(i)
                        residueSetE.add(j)
                        



    return residueSetA,residueSetE


parser = PDBParser(QUIET=True)
structure=parser.get_structure("6m0j","fixed_6m0j.pdb")

x=get_interface(structure,3)



