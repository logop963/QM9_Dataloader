
import os

import numpy as np

# Specify the path to the folder you want to iterate over
class FileReader:
    def __init__(self):
        self.smiles_string = []
        self.feature_values = []
    def smile_reader(self, folder_path):
        print("Reading SMILES...")
        for filename in os.listdir(folder_path):
            # Construct the full file path by joining the folder path and filename
            file_path = os.path.join(folder_path, filename)
            # Reading each file, getting each of the lines and extracting only letters for the SMILES string for rdkit
            with open(file_path, 'r') as file:
                smiles_line = file.readlines()[-2]
                smiles = smiles_line.split()[0]
                self.smiles_string.append(smiles)
        return self.smiles_string

    def feature_reader(self, folder_path):
        print("Reading features...")
        for filename in os.listdir(folder_path):
            # Construct the full file path by joining the folder path and filename
            file_path = os.path.join(folder_path, filename)
            # Reading each file, getting each of the lines and extracting only letters for the SMILES string for rdkit
            with open(file_path, 'r') as file:
                feature_line = file.readlines()[1]
                feature = feature_line.split()
                try:
                    feature_value = [float(feature[1+i]) for i in range(len(feature)-1)]
                except ValueError:
                    feature_value = None  # Handle non-numeric values
                self.feature_values.append(feature_value)
        return self.feature_values
