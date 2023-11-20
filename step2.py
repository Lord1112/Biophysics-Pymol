import step_2_modules
from Bio.PDB.NeighborSearch import NeighborSearch
import math
from Bio.PDB import PDBParser
from Bio.PDB import Structure

parser = PDBParser(QUIET=True)
st = parser.get_structure("cl_d6m0j.pdb", "6m0j.pdb")

# Assume you have two chains A and B
chain_A = st[0]["A"]
chain_B = st[0]["E"]
path_lib = 'parameters/aalib.lib'

res_lib = {}

with open(path_lib, 'r') as lib_file:
    for line in lib_file:
        # Skip lines starting with '#'
        if line.startswith('#'):
            continue
        # Split the line into columns
        columns = line.split()
        # Extract information from the columns
        residue = columns[0]
        atom = columns[1]
        atom_type = columns[2]
        charge = float(columns[3])  # Convert charge to float
        # Create a key for the dictionary (e.g., 'ALA_N')
        key = f'{residue}_{atom}'
        # Store the information in the dictionary
        res_lib[key] = {'residue': residue, 'atom': atom, 'atom_type': atom_type, 'charge': charge}

# Print or use the parsed data
txt_file_path = 'parameters/vdwprm.txt'

# Read lines from the file
with open(txt_file_path, 'r') as txt_file:
    ff_params = [line.strip() for line in txt_file]

# Print or use the parsed data
print()

l = step_2_modules.add_atom_parameters(st,res_lib,ff_params)

t = step_2_modules.calc_int_energies(chain_A, chain_B)





