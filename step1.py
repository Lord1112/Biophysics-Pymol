from Bio.PDB import PDBParser

def get_interface_residues(structure, chain, distance_threshold):
    interface_residues = set()

    for model in structure:
        for residue_a in model[chain]:
            for residue_e in model[chain]:
                if residue_a.resname != residue_e.resname and residue_a.get_vector() - residue_e < distance_threshold:
                    interface_residues.add(residue_a)
                    interface_residues.add(residue_e)

    return list(interface_residues)

# Cambia "tu_estructura.pdb" al nombre de tu archivo PDB y ajusta "A" y "B" a las cadenas relevantes.
parser = PDBParser(QUIET=True)
structure = parser.get_structure("6m0j", "hola.pdb")

#distancia de contacto que determinaste en PyMOL
distance_threshold = 2.0

# Obtenemos los residuos de la interfaz para cada cadena
interface_residues_chain_A = get_interface_residues(structure, "A", distance_threshold)
interface_residues_chain_E = get_interface_residues(structure, "E", distance_threshold)

# Imprimimos los resultados
print(f"Residuos de la interfaz en la cadena A: {interface_residues_chain_A}")
print(f"Residuos de la interfaz en la cadena E: {interface_residues_chain_E}")
