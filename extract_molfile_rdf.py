from rdkit import Chem
import rdkit.Chem.Lipinski as lip
import sys
import subprocess


class MolFileExtraction:

    end_regex = r'^M  END'
    molfile_no = r'^\$MFMT'

    start_index = []
    end_index = []

    def __init__(self, file: str) -> None:
        if file.endswith('.rdf'):
            self.file = file
            print(f'reading {self.file}')
            self.lines = open(self.file, 'r', encoding='utf-8').readlines()

            self.extract_indices()
        else:
            assert 'Wrong file type provided.'

    def extract_indices(self):

        if sys.platform == 'linux' or sys.platform == 'linux2' or sys.platform == 'darwin':
            start = subprocess.run(
                ['grep', '-n', f" '{self.molfile_no}' {self.file} > temp.temp"], capture_output=True)
            print(start.stdout.decode('utf-8').split('\r\n'))
            # added hint for developing on non-windows systems.

        elif sys.platform == 'win32':

            subprocess.run(
                ['powershell', '-Command', f"findstr -n '{self.molfile_no}' {self.file} > temp"], capture_output=True)

            self.mol_start = subprocess.run(
                ['powershell', '-Command', r"Get-Content .\temp | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

            # subprocess.run(['powershell', 'rm temp'])

            subprocess.run(
                ['powershell', '-Command', f"findstr -n '{self.end_regex}' {self.file} > temp2"], capture_output=True)

            self.mol_end = subprocess.run(
                ['powershell', '-Command', r"Get-Content .\temp2 | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

            self.molfiles = list(zip(self.mol_start, self.mol_end))

        self.calculate_HDHN()

    def calculate_HDHN(self):
        mol = ''.join(self.lines[int(self.molfiles[0][0]): int(self.molfiles[0][1])])

        m = Chem.rdmolfiles.MolFromMolBlock(mol)

        ND = lip.NumHDonors(m)
        NA = lip.NumHAcceptors(m)

        self.lines.insert(int(self.molfiles[0][1])+12, '$DTYPE MOL:RESULT(1):NumHDonors\n')
        self.lines.insert(int(self.molfiles[0][1])+13, f'$DATUM {ND}\n')
        self.lines.insert(int(self.molfiles[0][1])+14, '$DTYPE MOL:RESULT(1):NumHAcceptors\n')
        self.lines.insert(int(self.molfiles[0][1])+15, f'$DATUM {NA}\n')

        with open('mod_test_rdf.rdf', 'w') as f:
            f.writelines(''.join(self.lines))


obj = MolFileExtraction('test_rdf.rdf')
