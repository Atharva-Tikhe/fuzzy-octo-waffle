from rdkit import Chem
import rdkit.Chem.Lipinski as lip
import subprocess
import re
input_file = open('trial.rdf', 'r').readlines()

def process_pbm(filename):
    pbm_start = ' MOL:PBM_STRUCTURE$'
    pbm_end = ' MOL:PBM_FORMULA$'
    subprocess.run(
        ['powershell', '-Command', f"findstr -n -r '{pbm_start}' {filename} > ..\\advanced_testing\\temppbm"],
        capture_output=True)

    mol_start_pbm = subprocess.run(
        ['powershell', '-Command', r"Get-Content ..\advanced_testing\temppbm | %{$_.Split(':')[0]}"],
        capture_output=True).stdout.decode('utf-8').split("\r\n")

    subprocess.run(
        ['powershell', '-Command', f"findstr -n -r '{pbm_end}' {filename} > ..\\advanced_testing\\temppbm_end"],
        capture_output=True)

    mol_end_pbm = subprocess.run(
        ['powershell', '-Command', r"Get-Content ..\advanced_testing\temppbm_end | %{$_.Split(':')[0]}"],
        capture_output=True).stdout.decode('utf-8').split("\r\n")

    pbm_molfiles = list(zip(mol_start_pbm, mol_end_pbm))

    print(pbm_molfiles)
    for molfile in pbm_molfiles:
        if pbm_molfiles.index(molfile) == 0 and molfile[0] != '':
            mol = ''.join(input_file[int(molfile[0])+1:int(molfile[1])-1])

            m = Chem.rdmolfiles.MolFromMolBlock(mol)
            smiles = Chem.MolToSmiles(m, isomericSmiles=False)
            input_file.insert(int(molfile[1]) + 1, '$DTYPE MOL:PBM_SMILES\n')
            input_file.insert(int(molfile[1]) + 2, f'$DATUM {smiles}\n')
        elif molfile[0] != '':
            mol = ''.join(input_file[int(molfile[0]) + 3:int(molfile[1])+1])
            m = Chem.rdmolfiles.MolFromMolBlock(mol)
            smiles = Chem.MolToSmiles(m, isomericSmiles=False)
            input_file.insert(int(molfile[1]) + 1, '$DTYPE MOL:PBM_SMILES\n')
            input_file.insert(int(molfile[1]) + 2, f'$DATUM {smiles}\n')

    with open('output_rdf.rdf', 'w') as f:
        f.writelines(input_file)




process_pbm('trial.rdf')

'''
A
B
C
A
B
C
.
.
.
initial coords
Read all A and modify the file, save it
A updates coords
Read all B and modify the file, save it
B updates coords
Read all C and modify the file, save it
C updates coords (dont care)
exit

4 IOs
3 regex matches
index shifting 2,2,4
'''