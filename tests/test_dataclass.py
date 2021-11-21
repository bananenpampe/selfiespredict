from selfiespredict.data.load_data import Data_Cleaner #imports Data_cleaner
import unittest
from rdkit import Chem
import os
from copy import deepcopy
from rdkit import RDLogger

paper_smiles = "CNC(C)CC1=CC=C2C(=C1)OCO2"
paper_smiles_rdkit_canonical = "CNC(C)Cc1ccc2c(c1)OCO2"
paper_smiles_rdkit_canonical_tokenized = "C N C ( C ) C c 1 c c c 2 c ( c 1 ) O C O 2"
paper_smiles_tokenized = "C N C ( C ) C C 1 = C C = C 2 C ( = C 1 ) O C O 2"
paper_SELFIES = "[C][N][C][Branch1][C][C][C][C][=C][C][=C][C][=Branch1][Ring2][=C][Ring1][=Branch1][O][C][O][Ring1][=Branch1]"
paper_SELFIES_tokenized = "[C] [N] [C] [Branch1] [C][C][C][C][=C][C][=C][C][=Branch1][Ring2][=C][Ring1][=Branch1][O][C][O][Ring1][=Branch1]"

TEST_DATAPATH = "test_data/test_SMILES.txt"

TEST_DATAPATH = os.path.join(os.path.dirname(os.path.realpath(__file__)),TEST_DATAPATH)
RDLogger.DisableLog('rdApp.*')
#test if import if import wo path is None
#test if data is none after import
#test if raw_data

class TestDataClass(unittest.TestCase):
    def test_dataclass_none(self):
        """Checking if cleaner.raw_data and cleaner.data is None after empty init"""
        cleaner = Data_Cleaner()
        self.assertIsNone(cleaner.raw_data)
        self.assertIsNone(cleaner.data)

    def test_dataclass_afterload_none(self):
        """Checking if cleaner.data is None after init with PATH"""
        cleaner = Data_Cleaner(TEST_DATAPATH)
        self.assertIsNone(cleaner.data)
        self.assertIsInstance(cleaner.raw_data,list)

    def test_single_smile_raw_to_smile_tokenized(self):
        """Check if single smile gets correctly converted to SELFIE"""
        cleaner = Data_Cleaner()
        cleaner.raw_data = [paper_smiles, paper_smiles_rdkit_canonical]
        self.assertIsNone(cleaner.data)
        cleaner.data2smile()
        self.assertEqual(cleaner.data[0], paper_smiles_rdkit_canonical_tokenized)
        self.assertEqual(cleaner.data[1], paper_smiles_rdkit_canonical_tokenized)

    def test_single_smile_data_to_selfie(self):
        """Check if single smile (rdkit and SELFIES canonical)
        in .data gets correctly converted to SELFIE"""
        cleaner = Data_Cleaner()
        cleaner.data = [paper_smiles, paper_smiles_rdkit_canonical]
        cleaner.smile2selfie()
        self.assertEqual(cleaner.data[0], paper_SELFIES)
        self.assertEqual(cleaner.data[1], paper_SELFIES)

    def test_single_selfie_to_smile(self):
        """Check if single SELFIES gets tokenized correctly"""
        cleaner = Data_Cleaner()
        cleaner.data = [paper_SELFIES]
        cleaner.selfie2smile()
        self.assertEqual(cleaner.data[0], paper_smiles_rdkit_canonical_tokenized)

    def test_smiles_canonical_test_data(self):
        """Check SMILES on testdatset"""
        cleaner = Data_Cleaner(TEST_DATAPATH)
        self.assertTrue(all([repr1 == Chem.CanonSmiles(repr1) for repr1 in cleaner.raw_data]))

    def test_smiles_canonical_test_data(self):
        """Checks if smiles_tok -> selfies -> smiles_tok is rdkit canonical"""
        cleaner = Data_Cleaner(TEST_DATAPATH)
        cleaner.data2smile()
        this = deepcopy(cleaner.data)
        cleaner.smile2selfie()
        cleaner.selfie2smile()
        self.assertTrue(all([repr1[0] == repr1[1] for repr1 in zip(this,cleaner.data)]))




        #assertTrue(all([repr1 == conversion(repr1) for repr1 in data]))
