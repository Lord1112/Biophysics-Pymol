from Bio.PDB import PDBParser
from Bio.PDB import *

def get_interface(structure,threshold):
    chainA=structure[0]["A"]
    chainE=structure[0]["E"]
    residueSetA=set()
    residueSetE=set()
    for i in chainA:
        for y in i:
            for j in chainE:
                for z in j:                
                    if y-z <= threshold:
                        residueSetA.add(i)
                        residueSetE.add(j)
                        



    return residueSetA,residueSetE


parser = PDBParser(QUIET=True)
structure=parser.get_structure("6m0j","fixed_6m0j.pdb")

x=get_interface(structure,3)



