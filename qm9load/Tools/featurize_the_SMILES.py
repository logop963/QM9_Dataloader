import deepchem as dc
import pickle
import numpy as np

def rdkit_descriptors(smiles_data, make_hash=None):
    """
    Uses SMILES data and uses deepchem to featurize
    :param smiles_data: SMILES strings from X
    :return:
    """
    print("Calculating RDKit Fingerprint...")
    smiles_list = list(smiles_data)
    rdkit_vect = []
    rdkit_dict = {}
    for i, smiles in enumerate(smiles_list):
        featurizer = dc.feat.RDKitDescriptors()
        find_features = featurizer.featurize(str(smiles))
        rdkit_vect.append(find_features[0])
        if make_hash:
            rdkit_dict[smiles] = find_features[0]

    with open(
            r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Fingerprint_Dictionaries\RDKit_Fingerprint',
            'wb') as f:
        pickle.dump(rdkit_dict, f)
    return rdkit_vect


if __name__ == '__main__':
    with open(r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Dataset\QM9.pkl', 'rb') as f:
        data = pickle.load(f)
    smiles = data['SMILES String'].values
    output = rdkit_descriptors(smiles, True)
    print(output)