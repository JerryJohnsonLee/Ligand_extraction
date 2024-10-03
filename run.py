'''
Extract Ligands from a PDB File

Input - pdb_file: A file path to a PDB structure
peptide_ligand_length: Any peptide shorter than this specified length from the structure will also be treated as a ligand
remove_solvents_and_ions: By default set to True. If False, solvents and ions will also be returned
Output - A csv file containing the ligand names, their SMILES and the extracted sdf files

This tool extracts the ligand structures from a PDB file directly. 
It uses biopython to read in the structure and analyze the chain sequences. 
Any chain shorter than the specified peptide_ligand_length will be treated as a ligand. 
Otherwise, if a residue does not belong to standard amino acid code, it will be treated as a ligand 
    unless there is a standard amino acid in the same chain within 5 residues of the given position on both ends. 
Every residue is saved as a separate sdf file and the results are summarized into a data table.
'''

from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
import random
import string
import pandas as pd

# define standard_amino_acids
standard_amino_acids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
excluded_solvents_ions = ["HOH", "H2O"]
ligand_residue_check_window = 5

# Custom select class
class SelectProteinLigand(Select):
    def accept_residue(self, residue):
        # Check for the custom annotation in each residue
        if hasattr(residue, 'protein_and_ligand'):
            return residue.protein_and_ligand
        else:
            return True
        
def create_random_string(length=3):
    return ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))


def convert_pdb_to_sdf_and_smiles(identifier):
    mol = Chem.MolFromPDBFile(identifier + ".pdb")
    if mol is None:
        return ""
    else:
        Chem.MolToMolFile(mol, identifier + ".sdf")
        return Chem.MolToSmiles(mol)

def run(pdb_file, peptide_ligand_length=10, remove_solvents_and_ions=True):
    # read in the pdb file
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    ligands = []
    ligand_names = []


    # first save a structure without solvent and ions
    for residue in structure.get_residues():
        is_solvent_or_ion = False
        if residue.get_resname() in excluded_solvents_ions:
            is_solvent_or_ion = True
        else:
            # a residue that does not have any carbon atom is also considered as solvent or ion
            atoms_in_residue = [atom.element for atom in residue]
            if "C" not in atoms_in_residue:
                is_solvent_or_ion = True
        if is_solvent_or_ion and not remove_solvents_and_ions:
            ligands.append(residue)
            ligand_names.append(residue.get_resname())
        residue.protein_and_ligand = not is_solvent_or_ion
    
    # save proteins and ligands as a separate pdb file
    io = PDBIO()
    io.set_structure(structure)
    io.save("protein_and_ligand.pdb", select=SelectProteinLigand())

    # now further process the protein and ligand structure
    structure = parser.get_structure("structure", "protein_and_ligand.pdb")
    for model in structure:
        for chain in model:
            if len(chain) <= peptide_ligand_length:
                ligands.append(chain)
                ligand_names.append("peptide_"  + chain.get_id())
                continue
            else:
                for residue in chain:
                    if residue.get_resname() not in standard_amino_acids:
                        # check if there is a standard amino acid within 5 residues on both ends
                        check_before = False
                        for i in range(1, ligand_residue_check_window + 1):
                            try:
                                if chain[residue.id[1] - i].get_resname() in standard_amino_acids:
                                    check_before = True
                                    break
                            except:
                                continue

                        check_after = False
                        for i in range(1, ligand_residue_check_window + 1):
                            try:
                                if chain[residue.id[1] + i].get_resname() in standard_amino_acids:
                                    check_after = True
                                    break
                            except:
                                continue

                        if not (check_before and check_after):
                            ligands.append(residue)
                            ligand_names.append(residue.get_resname())

    # write out the ligands
    ligand_smiles = []
    ligand_files = []
    for ligand, ligand_name in zip(ligands, ligand_names):
        ligand_identifier = f"{ligand_name}_{create_random_string()}"
        io = PDBIO()
        io.set_structure(ligand)
        io.save(ligand_identifier + ".pdb")
        smiles = convert_pdb_to_sdf_and_smiles(ligand_identifier)
        ligand_files.append(ligand_identifier + ".sdf")
        ligand_smiles.append(smiles)

    # write out the results
    csv_filename = f"ligand_extraction_{create_random_string(length=6)}.csv"
    df = pd.DataFrame({"ligand_name": ligand_names, "SMILES": ligand_smiles, "SDF_path": ligand_files})
    # order according to SMILES length
    df = df.sort_values("SMILES", key=lambda x: x.str.len(), ascending=False)
    df.to_csv(csv_filename, index=False)

    print(f"The ligand information are reported in the file: {csv_filename}")
    print(f"The ligand sdf files are:", ligand_files)

if __name__ == "__main__":
    run("test_data/7rt4.pdb", remove_solvents_and_ions=False)