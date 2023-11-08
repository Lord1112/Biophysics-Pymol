from Bio.PDB import *
from Bio.PDB import Chain
from Bio.PDB import Model
from Bio.PDB import Structure
from Bio.PDB.PDBIO import PDBIO

def chainCleaner(chain,chainID):
    cleanedChain = Chain.Chain(chainID)
    for i in chain.get_residues():
        if i.resname in {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'}:
            cleanedChain.add(i)
    return(cleanedChain)

parser = PDBParser(QUIET=True)
unStructure=parser.get_structure("6m0j","6m0j.pdb")
unChainA = unStructure[0]["A"]
unChainE = unStructure[0]["E"]
        
clChainA= chainCleaner(unChainA,"A")
clChainE= chainCleaner(unChainE,"E")

clModel = Model.Model(0)

clModel.add(clChainA)
clModel.add(clChainE)

clStructure = Structure.Structure(0)

clStructure.add(clModel)

#Save the cleaned structure to a pdb file
pdbCreator=PDBIO()
pdbCreator.set_structure(clStructure)
pdbCreator.save("cl_d6m0j.pdb")

import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
base_dir_path=biobb_structure_checking.__path__[0]
base_path = ''
args = cts.set_defaults(base_dir_path,{'notebook':True})
args['input_structure_path'] = base_path + '6m0j.cif'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['debug'] = False
args['verbose'] = False
#Added missing defaults
args['time_limit'] = 3600
args['nocache'] = False
args['copy_input'] = None
args['build_warnings'] = False
args['coords_only'] = False
args['overwrite'] = False
args['output_format'] = 'pdb'
