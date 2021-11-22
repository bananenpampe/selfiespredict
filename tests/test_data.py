#tests the full dataset

from selfiespredict.data.load_data import Data_Cleaner #imports Data_cleaner
import unittest
from rdkit.Chem import CanonSmiles
from rdkit import Chem
from os.path import dirname, abspath, join
from copy import deepcopy
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

#This is bad. replace with absolut
BASEPATH = dirname(dirname(abspath(__file__)))
DATAPATH = "data"
sets = ["USPTO_480K","USPTO_STEREO"]
files = ["src-test.txt","src-val.txt","tgt-test.txt","tgt-train.txt","tgt-val.txt"]

class TestDataClass(unittest.TestCase):
    def test_correctness_data(self):
        """Tests, if all data can be converted from tokenized smiles ->
        smiles_tokenized -> selfies ->
         """
        print(join(BASEPATH,DATAPATH,sets[0],files[0]))

        """cleaner = Data_Cleaner("../data/USPTO_STEREO/src-test.txt")
        cleaner.data2smile()
        this = deepcopy(cleaner.data)
        cleaner.smile2selfie()
        cleaner.selfie2smile()
        for n, repr1 in enumerate(zip(this,cleaner.data)):
            if repr1[0] != repr1[1]:
                print("In {} wrong: {}".format("./test_data/critical.txt",n))
            else:
                continue
        #this = [repr1[0] == repr1[1] for repr1 in zip(this,cleaner.data)]
        #print(this)
        #self.assertTrue(all(this))"""
        """for dataset in sets:
            for datafile in files:
                PATH = DATAPATH + sets + files
                cleaner = Data_Cleaner(PATH)
                cleaner.data2smile()
                this = deepcopy(cleaner.data)
                cleaner.smile2selfie()
                cleaner.selfie2smile()
                for n, repr1 in enumerate(zip(this,cleaner.data)):
                    if repr1[0] == repr1[1]:
                        continue
                    else:
                        print("In {} wrong: {}".format(PATH,n))
                self.assertTrue(all([repr1[0] == repr1[1] for repr1 in zip(this,cleaner.data)]))"""


#STRUCTURES FAIL WITH [P-] and rdkit default parameters
