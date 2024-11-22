import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.decomposition import PCA
import numpy as np
import pickle
from utils import get_smiles_from_name

# Load the interaction data
with open('data/interaction_data_with_reduced_binding_features.pkl', 'rb') as file:
    interaction_data = pickle.load(file)

# PCA setup (on reduced fingerprints)
reduced_fingerprints = np.array([np.array(rf) for rf in interaction_data["reduced_fingerprint"]])
# print(f"Type of reduced_fingerprints: {type(reduced_fingerprints)}")
# print(f"Shape of reduced_fingerprints: {np.shape(reduced_fingerprints)}")
# print(f"Sample data from reduced_fingerprints: {reduced_fingerprints[:5]}")

pca = PCA(n_components=5)
pca.fit(reduced_fingerprints)

def query_openai_api(smell_description):
    """Query OpenAI to break down a smell description into a mix of compounds."""
    message = f"Break down the smell '{smell_description}' into a mix of 3-6 odorants."
    response = requests.post("https://api.openai.com/v1/completions", json={
        "model": "gpt-4", "prompt": message, "max_tokens": 100})
    if response.status_code == 200:
        return response.json()['choices'][0]['text']
    else:
        return {"error": f"Error querying OpenAI API: {response.text}"}

def generate_descriptors(smiles):
    """Generate molecular descriptors for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return np.array([
            Descriptors.MolWt(mol),
            Descriptors.TPSA(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumHDonors(mol)
        ])
    else:
        return None

def calculate_cartridge_activations(odorants):
    """Calculate cartridge activations for the odorant mix."""
    descriptors = {o['name']: generate_descriptors(get_smiles_from_name(o['name'])) for o in odorants}
    
    descriptor_projections = np.zeros(len(next(iter(descriptors.values()))))
    for odorant in odorants:
        descriptor = descriptors.get(odorant['name'])
        if descriptor is not None:
            weight = odorant['percentage'] / 100
            descriptor_projections += weight * descriptor

    cartridge_projections = normalize(descriptor_projections.reshape(1, -1), norm='l1') * 100
    return cartridge_projections[0].tolist()
