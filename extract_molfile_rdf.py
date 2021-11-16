from rdkit import Chem
import rdkit.Chem.Lipinski as lip
import re
import subprocess

class MolFileExtraction:

    start_regex = r'^  -ISIS- .*'
    end_regex = r'^M  END'
    molfile_no = r'^\$MFMT \$MIREG .*'


    start_index = []
    end_index = []

    def __init__(self, file: str) -> None:
        if file.endswith('.rdf'):
            self.lines = open('test.rdf', 'r', encoding='utf-8').readlines()
            # self.text = ''.join(open('test.rdf', 'r', encoding='utf-8').readlines())
            self.file = file
        else:
            assert 'Wrong file type provided.'

    def extract_indices(self):

        for line in self.lines:
            if line.startswith('$MFMT'):
                self.start_index.append(self.lines.index(line))
        check = re.compile(self.end_regex, re.MULTILINE)

        result = check.search(self.text)
        print(result)

        # if line.startswith('$DTYPE MOL:'):
        #     end_index.append(lines.index(line)-1)

        start = subprocess.run(['powershell', '-Command',"findstr -n '^\$MFMT $\n' test.rdf"], capture_output=True)
        print(start.stdout.decode('utf-8').split('\r\n'))

        molfiles = list(zip(self.start_index, self.end_index))
        print(molfiles)

        molfile_one = ''.join(self.lines[molfiles[0][0]+1: molfiles[0][1]+1])
        molfile_two = ''.join(self.lines[molfiles[1][0]+1: molfiles[1][1]+1])
        print(molfile_two)

        # with open('test.mol', 'w', encoding='utf-8') as f:
        #     f.writelines(molfile_one)
        # m = Chem.MolFromMolFile(r'path/to/*.mol')

        m = Chem.MolFromMolBlock(molfile_one)
        m1 = Chem.MolFromMolBlock(molfile_two)

        print(lip.NumHDonors(m))
        print(lip.NumHAcceptors(m))


        print(lip.NumHDonors(m1))
        print(lip.NumHAcceptors(m1))
