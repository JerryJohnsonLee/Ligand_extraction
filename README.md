Extract Ligands from a PDB File

Input - pdb_file: A file path to a PDB structure
peptide_ligand_length: Any peptide shorter than this specified length from the structure will also be treated as a ligand
remove_solvents_and_ions: By default set to True. If False, solvents and ions will also be returned
Output - A csv file containing the ligand names, their SMILES and the extracted sdf files

This tool extracts the ligand structures from a PDB file directly. It uses biopython to read in the structure and analyze the chain sequences. Any chain shorter than the specified peptide_ligand_length will be treated as a ligand. Otherwise, if a residue does not belong to standard amino acid code, it will be treated as a ligand unless there is a standard amino acid in the same chain within 5 residues of the given position on both ends. Every residue is saved as a separate sdf file and the results are summarized into a data table.

# Requirements
biopython
rdkit
pandas