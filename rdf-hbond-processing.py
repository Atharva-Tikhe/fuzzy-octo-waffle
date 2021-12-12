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
    results = []
    for molfile in pbm_molfiles:
        if molfile[0] != '':
            mol = ''.join(input_file[int(molfile[0])+1:int(molfile[1])-1])

            m = Chem.rdmolfiles.MolFromMolBlock(mol)
            smiles = Chem.MolToSmiles(m, isomericSmiles=False)
            results.append(int(molfile[1]))
            results.append(smiles)
    with open('output_rdf.rdf', 'w') as f:
        for line_index in range(len(input_file)):
            line = input_file[line_index]
            if line.endswith('MOL:PBM_FORMULA\n'):
                if line_index not in results:
                    line_index += 1
                input_file.insert(line_index + 1, f"$DTYPE MOL:PBM_SMILES\n{results[results.index(line_index) + 1]}\n")

            f.write(line)
        print(results)



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