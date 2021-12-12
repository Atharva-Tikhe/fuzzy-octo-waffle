from rdkit import Chem
import sys

molfile = sys.argv[1]

mol = Chem.rdmolfiles.MolFromMolBlock(molfile)

m = Chem.rdmolfiles.MolFromMolBlock(mol)

sys.stdout.write(f'{m}')
