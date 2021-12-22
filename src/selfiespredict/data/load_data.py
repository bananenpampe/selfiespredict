import selfies as sf
import selfiespredict
from selfiespredict.helpers.Helper_functions import smi_tokenizer
from rdkit.Chem import CanonSmiles
from rdkit import Chem
import os
import argparse
import gdown
import regex as re

def return_canonical(smiles):
    """Helper function that returns canolicalized SMILES with extended valency
    Depreciated as new RDKit version can handle valencies.
    """
    m = Chem.MolFromSmiles(smiles)
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

def selfie_splitter(input_data, brackets):
        """Helper function that splits (tokenizes) a list of SELFIES reaction smiles
        rule: "[Na+]" -> "[ Na + ]" if brackets
              "[Na+]" -> "Na +" if not brackets
        """
        output = [None]*len(input_data)
        for i in range(len(input_data)):
            filter_split = list(filter(None, re.split('(\W|\d)',input_data[i])))
            if not brackets:
                 filter_split = [item.replace(']','').replace('[','') for item in filter_split]
            removed_space = [i for i in filter_split if i != ""]

            if input_data[i] != "".join(removed_space) and brackets==False:
                  raise ValueError("ValueError exception thrown")

            string_list = ' '.join(removed_space)
            output[i] = string_list
        return output

class Data_Cleaner:
    """ Class for data preparation 'change made in colab2'

    Data_cleaner(FILEPATH)

    """

    def _load_data(self, DATAPATH):
        """
        Function that loads the dataset in .txt format from DATAPATH
        """
        with open(DATAPATH, 'r') as data_file:
            list_of_lists = []
            for line in data_file:
              stripped_line = line.strip()
              line_list = stripped_line.replace(" ", "")
              list_of_lists.append(CanonSmiles(line_list))

        self.raw_data = list_of_lists

    def __init__(self, DATAPATH=None, selfies_constraints=None):
        """Initializes Data_cleaner class.
           sets extended SELFIES constraints to deal with hypervalency
        """
        #adding selfies hypervalence constraints:

        if selfies_constraints is None:
            sf.set_semantic_constraints("hypervalent")
            constraints = sf.get_semantic_constraints()
            constraints['P-1'] = 7
            constraints['P'] = 6
            constraints['P+1'] = 5
            #14-18 should be max for organometallic transition metal complexes
            constraints['?'] = 18
            sf.set_semantic_constraints(constraints)
        else:
            sf.set_semantic_constraints(selfies_constraints)

        self.DATAPATH = DATAPATH
        self.raw_data = None
        self.data = None
        self.data_state = "raw"


        if self.DATAPATH:
            self._load_data(self.DATAPATH)

    def import_data(self,name):
        "Loads USPTO Stereo and USPTO480k into datafile as tokenized smiles from source"

        urls_fns_dict = {
        "USPTO_480k": [
            ("https://drive.google.com/uc?id=1RysNBvB2rsMP0Ap9XXi02XiiZkEXCrA8", "src-train.txt"),
            ("https://drive.google.com/uc?id=1CxxcVqtmOmHE2nhmqPFA6bilavzpcIlb", "tgt-train.txt"),
            ("https://drive.google.com/uc?id=1FFN1nz2yB4VwrpWaBuiBDzFzdX3ONBsy", "src-val.txt"),
            ("https://drive.google.com/uc?id=1pYCjWkYvgp1ZQ78EKQBArOvt_2P1KnmI", "tgt-val.txt"),
            ("https://drive.google.com/uc?id=10t6pHj9yR8Tp3kDvG0KMHl7Bt_TUbQ8W", "src-test.txt"),
            ("https://drive.google.com/uc?id=1FeGuiGuz0chVBRgePMu0pGJA4FVReA-b", "tgt-test.txt")
        ],
        "USPTO_STEREO": [
            ("https://drive.google.com/uc?id=1r3_7WMEor7-CgN34Foj-ET-uFco0fURU", "src-train.txt"),
            ("https://drive.google.com/uc?id=1HUBLDtqEQc6MQ-FZQqNhh2YBtdc63xdG", "tgt-train.txt"),
            ("https://drive.google.com/uc?id=1WwCH8ASgBM1yOmZe0cJ46bj6kPSYYIRc", "src-val.txt"),
            ("https://drive.google.com/uc?id=19OsSpXxWJ-XWuDwfG04VTYzcKAJ28MTw", "tgt-val.txt"),
            ("https://drive.google.com/uc?id=1FcbWZnyixhptaO6DIVjCjm_CeTomiCQJ", "src-test.txt"),
            ("https://drive.google.com/uc?id=1rVWvbmoVC90jyGml_t-r3NhaoWVVSKLe", "tgt-test.txt")
        ]
        }

        BASEPATH = os.path.dirname(os.path.dirname(os.path.dirname(selfiespredict.__file__)))
        data_path = os.path.join(BASEPATH,"data", "raw_data", name)

        os.makedirs(data_path, exist_ok=True)

        for url, fn in urls_fns_dict[name]:
            ofn = os.path.join(data_path, fn)
            if not os.path.exists(ofn):
                gdown.download(url, ofn, quiet=False)
                assert os.path.exists(ofn)
            else:
                print(f"{ofn} exists, skip downloading")

    def gen_txt(self,name):
        """
        Function generates txt file of tokenized SMILES and SELFIES for train/test/val
        """
        #generate tokenized selfies
        BASEPATH = os.path.dirname(os.path.dirname(os.path.dirname(selfiespredict.__file__)))
        data_path = os.path.join(BASEPATH,"data/tokenized_data/SELFIEwithBrackets", name)
        os.makedirs(data_path, exist_ok=True)
        raw_data_path = os.path.join("./data/raw_data", name)
        text_files = [f for f in os.listdir(raw_data_path) if f.endswith('.txt')]
        for file in text_files:
            data_file = os.path.join(raw_data_path, file)
            self._load_data(data_file)
            self.data2selfie()
            self.data = [' '.join(list(sf.split_selfies(item))) for item in self.data]
            output_file = os.path.join(data_path, file)
            with open(output_file, "w") as output:
                 for line in self.data:
                    output.write('%s\n' % line)

        #generate tokenized smiles
        data_path = os.path.join("./data/tokenized_data/SMILES", name)
        os.makedirs(data_path, exist_ok=True)
        raw_data_path = os.path.join("./data/raw_data", name)
        text_files = [f for f in os.listdir(raw_data_path) if f.endswith('.txt')]
        for file in text_files:
            data_file = os.path.join(raw_data_path, file)
            self._load_data(data_file)
            self.data2smile()
            output_file = os.path.join(data_path, file)
            with open(output_file, "w") as output:
                 for line in self.data:
                    output.write('%s\n' % line)

    def gen_SMILE_tokenized_SELFIES(self,name, brackets):
        """
        Function generates txt file of tokenized SMILES and SELFIES for train/test/val
        """
        #generate tokenized selfies
        BASEPATH = os.path.dirname(os.path.dirname(os.path.dirname(selfiespredict.__file__)))
        if brackets == True:
            data_path = os.path.join(BASEPATH,"data/tokenized_data/SMILE_tokenized_SELFIES_withBrackets", name)
        if brackets == False:
            data_path = os.path.join(BASEPATH,"data/tokenized_data/SMILE_tokenized_SELFIES_noBrackets", name)
        os.makedirs(data_path, exist_ok=True)
        raw_data_path = os.path.join("./data/raw_data", name)
        text_files = [f for f in os.listdir(raw_data_path) if f.endswith('.txt')]
        for file in text_files:
            data_file = os.path.join(raw_data_path, file)
            self._load_data(data_file)
            self.data2selfie()
            self.data = selfie_splitter(self.data, brackets)
            output_file = os.path.join(data_path, file)
            with open(output_file, "w") as output:
                 for line in self.data:
                    output.write('%s\n' % line)


    def data2smile(self):
        """
        Function that converts list entries to SMILE format
        """
        self.data = [smi_tokenizer(CanonSmiles(item)) for item in self.raw_data]
        self.data_state = "tokenized SMILE"

    def data2selfie(self):
        """
        Function that converts list entries to SELFIE format
        """
        self.data = [sf.encoder(item,strict=False) for item in self.raw_data]
        self.data_state = "SELFIE"

    def smile2selfie(self):
        """
        Function that converts SMILE list entries to SELFIE format
        """
        if self.data_state == "SMILE":
            self.data = [sf.encoder(item.replace(" ", ""),strict=False) for item in self.data]
            self.data_state = "SELFIE"
        else:
            print("Data in wrong state")

    def selfie2smile(self):
        """
        Function that converts SELFIE list entries to SMILE format
        """
        if self.data_state == "SELFIE":
            self.data = [smi_tokenizer(CanonSmiles(sf.decoder(item))) for item in self.data]
            self.data_state = "SMILE"
        else:
            print("Data in wrong state")

    def tokenize_selfie(self):
        """Function that tokenizes SELFIES
           Rule: [Na+][C-] -> Na+ C-
        """
        if self.data_state == "SELFIE":
            self.data = [' '.join(list(sf.split_selfies(item))).replace("[", "").replace("]", "") for item in self.data]
            self.data_state = "tokenized SELFIE"
        else:
            print("Data in wrong state")
