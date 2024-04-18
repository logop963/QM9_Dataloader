import deepchem as dc
import pickle
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def morgan_featurizer(smiles_data, make_hash=None):
    """
    Uses SMILES data and uses deepchem to featurize
    :param smiles_data: SMILES strings from X
    :return:
    """
    print("Calculating Morgan Fingerprint...")
    smiles_list = list(smiles_data)
    morgan_vect = []
    morgan_dict = {}
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        morgan = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2)
        morgan_str = morgan.ToBitString()
        morgan_array = np.array([int(bit) for bit in morgan_str])
        print(f'Morgan Shape: {morgan_array.shape}')
        morgan_vect.append(morgan_array)
        if make_hash:

            morgan_dict[smiles] = morgan_array

    with open(r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Fingerprint_Dictionaries\Morgan_Vector',
              'wb') as f:
        pickle.dump(morgan_dict, f)
    return morgan_vect


if __name__ == "__main__":
    with open(r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Dataset\QM9.pkl', 'rb') as f:
        data = pickle.load(f)
    smiles = data['SMILES String'].values
    output = morgan_featurizer(smiles, True)

    print(output)