from rdkit import Chem
from Bio import PDB

def run(input_pdb: str):
  # Load the PDB file
  parser = PDB.PDBParser()
  structure = parser.get_structure('protein', input_pdb)

  # Check if the protein has a single model and single chain
  if len(structure) != 1:
    raise ValueError('The protein has more than one model')
  if len(list(structure[0].get_chains())) != 1:
    raise ValueError('The protein has more than one chain')

  # Convert PDB to SDF
  mol = Chem.rdmolfiles.MolFromPDBFile(input_pdb)
  sdf_filename = input_pdb.replace('.pdb', '.sdf')
  Chem.rdmolfiles.MolToMolFile(mol, sdf_filename)

  return sdf_filename