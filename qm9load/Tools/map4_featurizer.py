import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

def maccs_featurizer(smiles_data, make_hash=None):
    """
    Uses SMILES data and MACCS fingerprinting to generate fingerprints
    :param smiles_data: SMILES strings from X
    :return: MACCS fingerprints as numpy array
    """
    print("Calculating MAACS Fingerprint...")
    smiles_list = list(smiles_data)
    dim = 166  # Dimensions of the MACCS fingerprint
    fingerprint_list = []
    maccs_dict = {}
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Generate MACCS keys fingerprint directly as a bit vector
            fp = AllChem.GetMACCSKeysFingerprint(mol)
            arr = np.array(fp, dtype=int)
            fingerprint_list.append(arr)
        else:
            # Handle invalid SMILES
            print("Invalid SMILES:", smiles)
            fingerprint_list.append(np.zeros(dim, dtype=int))  # Use zero vector for invalid SMILES
        if make_hash:
            maccs_dict[smiles] = arr
    with open(
            r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Fingerprint_Dictionaries\MACCS_Fingerprint',
            'wb') as f:
        pickle.dump(maccs_dict, f)

    return fingerprint_list

if __name__ == "__main__":
    with open(r'C:\Users\logop\PycharmProjects\MolecularAIProject\Peptide-ML\Dataset\QM9.pkl', 'rb') as f:
        data = pickle.load(f)
    smiles = data['SMILES String'].values
    output = maccs_featurizer(smiles, True)
    print(output)
