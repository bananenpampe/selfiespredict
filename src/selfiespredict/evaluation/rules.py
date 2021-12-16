from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
import selfies as sf
import numpy as np

def bond_count_difference_rule(mols):
    return mols[1].GetNumBonds() - mols[0].GetNumBonds()

def do_nothing_rule(mols):
    return 0

def return_rule_mat(rules,n_samples):
    return np.zeros((n_samples,len(rules)+2))


def return_feature_mat(PATH_EDUCT,PATH_PRODUCT,PATH_TRANS,rules=[],repr_type="SMILES"):

    sf.set_semantic_constraints("hypervalent")
    constraints = sf.get_semantic_constraints()
    constraints['P-1'] = 7
    constraints['P'] = 6
    constraints['P+1'] = 5
    #14-18 should be max for organometallic transition metal complexes
    constraints['?'] = 18
    sf.set_semantic_constraints(constraints)

    with open(PATH_EDUCT,"r") as e, open(PATH_PRODUCT,"r") as p, open(PATH_TRANS,"r") as t:
        lengths = [sum(1 for line in file) for file in [e,p,t]]
        if not lengths[0] == lengths[1] == lengths[2]:
            raise AssertionError("Something went wrong")

    with open(PATH_EDUCT,"r") as e, open(PATH_PRODUCT,"r") as p, open(PATH_TRANS,"r") as t:
        rule_matrix = return_rule_mat(rules,lengths[0])

        for n, lines in enumerate(zip(e,p,t)):
            #print("yo")
            lines = [ line.replace(' ', '').rstrip("\n") for line in lines]

            if repr_type == "SMILES":

                representations = []

                for smiles in lines:
                    try:
                        representations.append(Chem.MolFromSmiles(smiles))
                    except:
                        representations.append("NOTVALID")
                        print("smiles not valid")

            elif repr_type == "SELFIES":

                representations = []

                for selfies_repr in lines[:2]:
                    try:
                        representations.append(Chem.MolFromSmiles(sf.decoder(selfies_repr)))
                    except:
                        representations.append("NOTVALID")

            rule_matrix[n][0] = n

            if lines[2] == lines[1]:
                rule_matrix[n][1] = 1
            else:
                rule_matrix[n][1] = 0

            for j,rule in enumerate(rules):
                rule_matrix[n][j+2] = rule(representations)

    return rule_matrix

if __name__ == "__main__":
    THIS = return_feature_mat("./data/tokenized_data/SELFIE/USPTO_480k/src-val.txt","./data/tokenized_data/SELFIE/USPTO_480k/tgt-val.txt","./data/tokenized_data/SELFIE/USPTO_480k/predictions_best1_selfies_with_brackets.txt",rules=[bond_count_difference_rule],repr_type="SELFIES")
    #prints ERROR
    print(np.sum(THIS[:,1])/len(THIS))
