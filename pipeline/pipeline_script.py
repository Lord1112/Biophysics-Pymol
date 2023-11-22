from Bio.PDB import *
from Bio.PDB import Chain
from Bio.PDB import Model
from Bio.PDB import Structure
from Bio.PDB.PDBIO import PDBIO
import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
import argparse, sys, os, math
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select

def chainCleaner(chain,chainID):
    cleanedChain = Chain.Chain(chainID)
    for i in chain.get_residues():
        if i.resname in {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'}:
            cleanedChain.add(i)
    return(cleanedChain)

cleaningPath=sys.argv[1]

parser = PDBParser(get_header=True,QUIET=True)
structure=parser.get_structure("6m0j",cleaningPath)


unChainA = structure[0]["A"]
unChainE = structure[0]["B"]

clChainA= chainCleaner(unChainA,"A")
clChainE= chainCleaner(unChainE,"B")

pdbCreator=PDBIO()
clModel=Model.Model(0)
clStructure=Structure.Structure(0)
clModel.add(clChainA)
clModel.add(clChainE)
clStructure.add(clModel)
pdbCreator.set_structure(clStructure)
pdbCreator.save("clean.pdb")


base_dir_path=biobb_structure_checking.__path__[0]
base_path = ''

args = cts.set_defaults(base_dir_path,{'notebook':True})
args['input_structure_path'] = base_path + '6m0j.cif'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['input_structure_path'] = 'clean.pdb'
args['output_structure_path'] = 'cl_dm0j_fixed.pdb'
args['output_structure_path_charges'] = 'cl_d6m0j_fixed.pdbqt'
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

st_c = StructureChecking(base_dir_path, args)

st_c.chiral()
st_c.fixside()
st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')
st_c.add_hydrogen('auto')
st_c._save_structure("fixed.pdb")



# Functions needed to import the parameters for VanderWaals or the residue library
class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()

# Loading residue library from aaLib.lib
residue_library = ResiduesDataLib('parameters/aalib.lib')

# Loading VdW parameters
ff_params = VdwParamset('parameters/vdwprm.txt')

# Set the pdb_path and load the structure
pdb_path = "fixed.pdb"
# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)
# Loading structure
st = parser.get_structure('STR', pdb_path)

# Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

# Calculate electrostatic interactions
def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

# Calculate Van der Waals interactions
def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

# Calculate the solvation
def calc_solvation(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala

# Add atom parameters
def add_atom_parameters(st, res_lib, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

add_atom_parameters(st, residue_library, ff_params)

# Get the interface residues which are below a given distance.
def get_interface(st, dist):
    ''' Detects interface residues within a distance(dist)
        Assumes two chains, i.e. a unique interface set per chain.
    '''
    select_ats = []
    for at in st.get_atoms():
        # Skip Hydrogens to reduce time
        if at.element != 'H':
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {}
    # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        #Only different chains
        res1 = at1.get_parent()
        ch1 = res1.get_parent()
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface

get_interface(st, 3.5)

def calc_int_energies(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala

add_atom_parameters(st, residue_library, ff_params)
#COmplete path to nacess
NACCESS_BINARY = os.getcwd()+"/soft/NACCESS/naccess"
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

io = PDBIO()
st_chains = {}

# Get the chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0
        
for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

# Set a maximum distance (cut-off distance)
max_dist = 3.5

# Distance bigger than 0 gets the residues on interface.
if max_dist > 0.:
    interface = get_interface(st, max_dist)

# Inizialize dictionaries to store energies individually for each chain
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

# Inizialize conter variables to 0
totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {}

# Get the chains IDs, which won't be always A or B.
# Inizialize a variable for the solvation of each chain
chids = []
for ch in st[0]:
    chids.append(ch.id)
    totalSolvMon[ch.id] = 0

# Total counter for the binding energy
total = 0.

# Calculate the interaction energies for residues between chains A and E
# The output will be written in the int_energy_DG file

for ch in st[0]:
    for res in ch.get_residues():
        if max_dist > 0 and res not in interface[ch.id]:
            continue

        # Get the energies
        elec[res], elec_ala[res], vdw[res], vdw_ala[res] = calc_int_energies(st[0], res)
        solvAB[res], solvAB_ala[res] = calc_solvation(st[0], res)
        solvA[res], solvA_ala[res] = calc_solvation(st_chains[ch.id], st_chains[ch.id][0][ch.id][res.id[1]])

        # Add the values
        totalIntElec += elec[res]
        totalIntVdw += vdw[res]
        totalSolv += solvAB[res]
        totalSolvMon[ch.id] += solvA[res]

        total += elec[res] + vdw[res] + solvAB[res] - solvA[res]

print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
print('{:20}: {:11.4f}'.format('DG int AB-A-B', total))
