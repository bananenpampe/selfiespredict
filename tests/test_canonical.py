import unittest
from rdkit import Chem

def return_canonical(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

class TestCanonicalDemo(unittest.TestCase):
    def test_equivalent(self):
        """Is c1ccccc1 rdkit canonical"""
        benzene_smiles_can = "c1ccccc1"
        converted = return_canonical(benzene_smiles_can)
        self.assertEqual(benzene_smiles_can, converted, "should be canonical")

    @unittest.expectedFailure
    def test_NON_equivalent(self):
        """Is C1=CC=CC=C1 rdkit canonical"""
        benzene_smiles_other = "C1=CC=CC=C1"
        converted = return_canonical(benzene_smiles_other)
        self.assertEqual(benzene_smiles_other, converted)






if __name__ == '__main__':
    unittest.main()
