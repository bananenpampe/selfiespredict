#tests the full dataset

from selfiespredict.data.load_data import Data_Cleaner #imports Data_cleaner
import unittest
from rdkit.Chem import CanonSmiles
from rdkit import Chem
import os
from copy import deepcopy
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

DATAPATH = "../data/"
sets = ["USPTO_480K/","USPTO_STEREO/"]
files = ["src-test.txt","src-val.txt","tgt-test.txt","tgt-train.txt","tgt-val.txt"]

class TestDataClass(unittest.TestCase):
    def test_correctness_data(self):
        """Tests, if all data can be converted from tokenized smiles ->
        smiles_tokenized -> selfies ->
         """
        cleaner = Data_Cleaner("../data/USPTO_480K/src-test.txt")
        cleaner.data2smile()
        this = deepcopy(cleaner.data)
        cleaner.smile2selfie()
        cleaner.selfie2smile()

        for n, repr1 in enumerate(zip(this,cleaner.data)):
            if repr1[0] == repr1[1]:
                continue
            else:
                print(n)

        self.assertTrue(all([repr1[0] == repr1[1] for repr1 in zip(this,cleaner.raw_data)]))

#STRUCTURES FAIL WITH [P-] and rdkit default parameters
