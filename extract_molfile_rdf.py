from rdkit import Chem
import rdkit.Chem.Lipinski as lip
import sys
import subprocess


class MolFileExtraction:

    end_regex = '  END$'
    end_regex_protac = 'MOL:PD_FORMULA$'
    # molfile_no = r'^\$MFMT'
    molfile_no = '$DTYPE MOL:PD_STRUCTURE'

    UBM_start = ' MOL:UBM_STRUCTURE$'
    UBM_end = ' MOL:UBM_FORMULA$'

    PBM_start = ' MOL:PBM_STRUCTURE$'
    PBM_end = ' MOL:PBM_FORMULA$'
    start_index = []
    end_index = []

    def __init__(self, file: str) -> None:
        if file.endswith('.rdf'):
            self.file = file
            print(f'reading {self.file}')
            self.lines = open(self.file, 'r', encoding='utf-8').readlines()
            self.insert_list = []
            self.insert_list_pbm = []
            self.insert_list_ubm = []
            # self.extract_indices()
            self.insert_smiles_pbm()
            self.insert_smiles_ubm()
            self.extract_protac()

        else:
            assert 'Wrong file type provided.'

    def insert_smiles_pbm(self):
        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.PBM_start}' {self.file} > advanced_testing\\temppbm"], capture_output=True)

        self.mol_start_pbm = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\advanced_testing\temppbm | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.PBM_end}' {self.file} > advanced_testing\\temppbm_end"], capture_output=True)

        self.mol_end_pbm = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\advanced_testing\temppbm_end | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        self.molfiles_pbm = list(zip(self.mol_start_pbm, self.mol_end_pbm))

        for molblock in range(len(self.molfiles_pbm)):
            if self.molfiles_pbm[molblock][0] != '':
                mol = ''.join(self.lines[int(
                    self.molfiles_pbm[molblock][0])+1: int(self.molfiles_pbm[molblock][1])-1])
                print(mol)
                m = Chem.rdmolfiles.MolFromMolBlock(mol)

                SMILES = Chem.MolToSmiles(m, isomericSmiles=False)

                insert = [self.molfiles_pbm[molblock][1], SMILES]
                self.insert_list_pbm.append(insert)
                print(self.insert_list_pbm)

        for insert_info in self.insert_list_pbm:
            if self.insert_list_pbm.index(insert_info) == 0:
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:UBM_SMILES\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')

            else:
                insert_info[0] = int(insert_info[0]) + 2
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:UBM_SMILES\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')

    def insert_smiles_ubm(self):
        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.UBM_start}' {self.file} > advanced_testing\\tempubm"], capture_output=True)

        self.mol_start_ubm = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\advanced_testing\tempubm | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.UBM_end}' {self.file} > advanced_testing\\tempubm_end"], capture_output=True)

        self.mol_end_ubm = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\advanced_testing\tempubm_end | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        self.molfiles_ubm = list(zip(self.mol_start_ubm, self.mol_end_ubm))

        for molblock in range(len(self.molfiles_ubm)):
            if self.molfiles_ubm[molblock][0] != '':
                mol = ''.join(self.lines[int(
                    self.molfiles_ubm[molblock][0])+3: int(self.molfiles_ubm[molblock][1])+1])
                print(mol)

                m = Chem.rdmolfiles.MolFromMolBlock(mol)

                SMILES = Chem.MolToSmiles(m, isomericSmiles=False)

                insert = [self.molfiles_ubm[molblock][1], SMILES]
                self.insert_list_ubm.append(insert)


        for insert_info in self.insert_list_ubm:
            if self.insert_list.index(insert_info) == 0:
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:UBM_SMILES\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')

            else:
                insert_info[0] = int(insert_info[0]) + 4
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:UBM_SMILES\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')

    def extract_protac(self):
        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.molfile_no}' {self.file} > temp"], capture_output=True)

        self.mol_start = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\temp | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        subprocess.run(
            ['powershell', '-Command', f"findstr -n -r '{self.end_regex_protac}' {self.file} > temp2"], capture_output=True)

        self.mol_end = subprocess.run(
            ['powershell', '-Command', r"Get-Content .\temp2 | %{$_.Split(':')[0]}"], capture_output=True).stdout.decode('utf-8').split("\r\n")

        self.molfiles = list(zip(self.mol_start, self.mol_end))

        for molblock in range(len(self.molfiles)):
            if self.molfiles[molblock][0] != '':

                mol = ''.join(self.lines[int(
                    self.molfiles[molblock][0])+6: int(self.molfiles[molblock][1])-6])

                print(mol)

                m = Chem.rdmolfiles.MolFromMolBlock(mol)

                ND = lip.NumHDonors(m)
                NA = lip.NumHAcceptors(m)
                SMILES = Chem.MolToSmiles(m, isomericSmiles=False)

                insert = [self.molfiles[molblock][1], ND, NA, SMILES]
                self.insert_list.append(insert)

        for insert_info in self.insert_list:
            if self.insert_list.index(insert_info) == 0:
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:PD_HYDROGEN_BOND_DONOR\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')
                self.lines.insert(
                    int(insert_info[0])+3, '$DTYPE MOL:PD_HYDROGEN_BOND_ACCEPTOR\n')
                self.lines.insert(
                    int(insert_info[0])+4, f'$DATUM {insert_info[2]}\n')
                self.lines.insert(
                    int(insert_info[0])+5, f'$DTYPE MOL:PD_SMILES\n')

                self.lines.insert(
                    int(insert_info[0])+6, f'$DATUM {insert_info[3]}\n')

            else:
                insert_info[0] = int(insert_info[0]) + 6
                print('processing second protac')
                print(insert_info)
                self.lines.insert(
                    int(insert_info[0])+1, '$DTYPE MOL:PD_HYDROGEN_BOND_DONOR\n')
                self.lines.insert(
                    int(insert_info[0])+2, f'$DATUM {insert_info[1]}\n')
                self.lines.insert(
                    int(insert_info[0])+3, '$DTYPE MOL:PD_HYDROGEN_BOND_ACCEPTOR\n')
                self.lines.insert(
                    int(insert_info[0])+4, f'$DATUM {insert_info[2]}\n')
                self.lines.insert(
                    int(insert_info[0])+5, f'$DTYPE MOL:PD_SMILES\n')
                self.lines.insert(
                    int(insert_info[0])+6, f'$DATUM {insert_info[3]}\n')

        print('writing file')
        with open(r'advanced_testing\mod_test_rdf.rdf', 'w') as f:
            f.writelines(''.join(self.lines))


obj = MolFileExtraction(r'advanced_testing\trial.rdf')
