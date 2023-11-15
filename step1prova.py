from Bio.PDB import PDBParser

def calculate_distance(atom1, atom2):
    # Calculate the Euclidean distance between two atoms
    return atom1 - atom2

def find_interface_residues(chain1, chain2, distance_threshold):
    interface_residues_chain1 = set()
    interface_residues_chain2 = set()

    for residue1 in chain1:
        for atom1 in residue1:
            for residue2 in chain2:
                for atom2 in residue2:
                    distance = calculate_distance(atom1, atom2)
                    if distance < distance_threshold:
                        # Store a tuple with residue index and residue object
                        interface_residues_chain1.add((residue1.id[1], residue1))
                        interface_residues_chain2.add((residue2.id[1], residue2))

    return interface_residues_chain1, interface_residues_chain2

# Load the cleaned structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("cleaned_structure", "6m0j.pdb")

# Assume you have two chains A and B
chain_A = structure[0]["A"]
chain_B = structure[0]["E"]

# Define the distance threshold for interface residues
distance_threshold = 3.0  

# Find  thie interface residues
interface_residues_A, interface_residues_B = find_interface_residues(chain_A, chain_B, distance_threshold)

print("Interface residues in Chain A:")
for residue_index, residue_obj in interface_residues_A:
    print(f"Residue {residue_index} in Chain A: {residue_obj}")

print("\nInterface residues in Chain B:")
for residue_index, residue_obj in interface_residues_B:
    print(f"Residue {residue_index} in Chain B: {residue_obj}")

