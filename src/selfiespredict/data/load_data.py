import selfies as sf
from selfiespredict.helpers.Helper_functions import smi_tokenizer
from rdkit.Chem import CanonSmiles
from rdkit import Chem

def return_canonical(smiles):
    m = Chem.MolFromSmiles(smiles)
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m,Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_KEKULIZE|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


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

    def __init__(self, DATAPATH=None):

        self.DATAPATH = DATAPATH
        self.raw_data = None
        self.data = None


        if self.DATAPATH:
            self._load_data(self.DATAPATH)

    def data2smile(self):
        """
        Function that converts list entries to SMILE format
        """
        self.data = [smi_tokenizer(return_canonical(item)) for item in self.raw_data]

    def data2selfie(self):
        """
        Function that converts list entries to SELFIE format
        """
        self.data = [sf.encoder(item,strict=False) for item in self.raw_data]

    def smile2selfie(self):
        """
        Function that converts SMILE list entries to SELFIE format
        """

        self.data = [sf.encoder(item.replace(" ", ""),strict=False) for item in self.data]

    def selfie2smile(self):
        """
        Function that converts SELFIE list entries to SMILE format
        """

        self.data = [smi_tokenizer(return_canonical(sf.decoder(item))) for item in self.data]
