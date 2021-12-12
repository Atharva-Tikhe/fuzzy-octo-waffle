from rdkit import Chem
import rdkit.Chem.Lipinski as lip
import sys

molfile = sys.argv[1]

mol = Chem.rdmolfiles.MolFromMolBlock(molfile)

ND = lip.NumHDonors(mol)
NA = lip.NumHAcceptors(mol)

sys.stdout.write(f'{NA}|{ND}')
