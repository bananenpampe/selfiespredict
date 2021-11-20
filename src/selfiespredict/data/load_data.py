import selfies as sf
from selfiespredict.helpers.Helper_functions import smi_tokenizer    

class Data_Cleaner:
    """ Class for data preparation 'change made in colab2'

    Data_cleaner(FILEPATH)

    """

    def _load_data(self, DATAPATH):
        """
        Function that loads the dataset in .txt format from DATAPATH
        """
        data_file = open(DATAPATH, 'r')

        list_of_lists = []
        for line in data_file:
          stripped_line = line.strip()
          line_list = stripped_line.replace(" ", "")
          list_of_lists.append(line_list)

        self.raw_data = list_of_lists

    def __init__(self, DATAPATH=None):

        self.DATAPATH = DATAPATH
        self.raw_data = None
        self.data = None
        sf.set_semantic_constraints("hypervalent")

        if self.DATAPATH:
            self._load_data(self.DATAPATH)

    def data2smile(self):
        """
        Function that converts list entries to SMILE format
        """
        self.data = [smi_tokenizer(item) for item in self.raw_data]

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
        self.data = [smi_tokenizer(sf.decoder(item)) for item in self.data]
