import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.decomposition import PCA
import numpy as np
import pickle
from utils import get_smiles_from_name
from sklearn.preprocessing import normalize
from config import OPENAI_API_KEY
import os
from openai import OpenAI, APIConnectionError, RateLimitError, BadRequestError, APIError
import json
import re


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

# Initialize OpenAI client
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

def query_openai_api(smell_description):
    """
    Queries the OpenAI API with a smell description to retrieve a compound mix.

    Parameters:
        smell_description (str): Description of the smell to analyze.

    Returns:
        dict: JSON object containing the compound mix.
    """
    try:
        # Construct the message for the API call
        # message = (
        #     f"You are an expert chemist. Always respond with concise JSON data."
        #     f"Break down the smell '{smell_description}' into a mix of 3-6 odorants "
        #     f"with percentages summing to 100%. Only output valid JSON without any explanation or additional text: "
        #     f'{{ "odorants": [ {{"name": "odorant1", "percentage": 40}}, ... ] }}'
        # )
        message = (
            f"Break down the smell '{smell_description}' into a mix of 3-6 odorants (specific chemical compounds, if varied choose a specific compound) with percentages summing to 100%. Format the output as JSON: {{ 'odorants': [ {{'name': 'odorant1', 'percentage': 40}}, ... ] }}"
        )
        

        # Send the API request
        print("Sending request to OpenAI API...")
        response = client.chat.completions.create(
            model="gpt-4o-mini",  # Replace with the correct model
            messages=[{"role": "system", "content": message}],
            max_tokens=300,
            temperature=0.1,
        )

        # Extract and return the JSON content
        if response and hasattr(response, 'choices') and response.choices:
            print(response)
            return response.choices[0].message.content.strip()
        else:
            print("No valid 'choices' found in the response.")
            return None

    except APIConnectionError as e:
        print(f"API Connection Error: {e}")
        return {"error": "A connection error occurred with the OpenAI API."}
    except RateLimitError as e:
        print(f"Rate Limit Error: {e}")
        return {"error": "Rate limit exceeded. Please try again later."}
    except BadRequestError as e:
        print(f"Bad Request Error: {e}")
        return {"error": "Invalid request to OpenAI API. Check the inputs."}
    except APIError as e:
        print(f"OpenAI API Error: {e}")
        return {"error": "An error occurred with the OpenAI API."}
    except Exception as e:
        print(f"Unexpected Error: {e}")
        return {"error": f"Unexpected error: {str(e)}"}

def generate_descriptors(smiles):
    """Generate molecular descriptors for a given SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            descriptors = [
                Descriptors.MolWt(mol),            # Molecular weight
                Descriptors.TPSA(mol),             # Topological Polar Surface Area
                Descriptors.MolLogP(mol),          # LogP (octanol-water partition coefficient)
                Descriptors.NumHAcceptors(mol),    # Number of hydrogen bond acceptors
                Descriptors.NumHDonors(mol)        # Number of hydrogen bond donors
            ]
            return np.array(descriptors)
        else:
            print(f"Invalid SMILES: {smiles}")
            return None
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return None

# def calculate_cartridge_activations(odorants):
#     """
#     Process a list of odorants with their percentages to calculate cartridge activations.

#     Parameters:
#         odorants (list): List of dictionaries, each containing "name" and "percentage".

#     Returns:
#         dict: Cartridge activations and approximation accuracy.
#     """
#     # Step 1: Convert each compound to SMILES
#     smiles_list = {o['name']: get_smiles_from_name(o['name']) for o in odorants}
    
#     # Step 2: Generate descriptors for each compound
#     descriptors = {
#         name: generate_descriptors(smiles)
#         for name, smiles in smiles_list.items()
#         if smiles
#     }
    
#     if not descriptors:
#         raise ValueError("No valid SMILES or descriptors generated.")
    
#     # Step 3: Represent input as a weighted combination of primary cartridges
#     cartridge_projections = np.zeros(10)  # Initialize a 10-dimensional vector for the projection

#     for o in odorants:
#         name = o['name']
#         weight = o['percentage'] / 100
#         if name in descriptors and descriptors[name] is not None:
#             descriptor = descriptors[name][:5]  # Use only the first 5 descriptors
#             cartridge_projections += weight * descriptor

#     # Normalize cartridge projections
#     cartridge_percentages = normalize(cartridge_projections.reshape(1, -1), norm='l1') * 100

#     # Step 4: Calculate approximation accuracy
#     variance_explained = np.sum(pca.explained_variance_ratio_[:5]) * 100

#     return {
#         "cartridge_activations": cartridge_percentages[0].tolist(),
#         "approximation_accuracy": variance_explained
#     }

def calculate_cartridge_activations(odorants):
    """
    Process a list of odorants with their percentages to calculate cartridge activations.

    Parameters:
        odorants (list): List of dictionaries, each containing "name" and "percentage".

    Returns:
        dict: Cartridge activations and approximation accuracy.
    """
    # Step 1: Convert each compound to SMILES
    smiles_list = {o['name']: get_smiles_from_name(o['name']) for o in odorants}
    
    # Step 2: Generate descriptors for each compound
    descriptors = {
        name: generate_descriptors(smiles)
        for name, smiles in smiles_list.items()
        if smiles
    }
    
    if not descriptors:
        raise ValueError("No valid SMILES or descriptors generated.")
    
    # Step 3: Represent input as a weighted combination of primary cartridges
    descriptor_projections = np.zeros(pca.n_features_in_)  # Match PCA dimensionality

    for o in odorants:
        name = o['name']
        weight = o['percentage'] / 100
        if name in descriptors and descriptors[name] is not None:
            descriptor = descriptors[name]
            descriptor_projections[:len(descriptor)] += weight * descriptor

    # Step 4: Reduce dimensions using PCA
    reduced_descriptor_projections = pca.transform(descriptor_projections.reshape(1, -1))[0]

    # Step 5: Split into positive and negative projections
    positive_projections = np.zeros(5)
    negative_projections = np.zeros(5)

    for i in range(5):
        if reduced_descriptor_projections[i] > 0:
            positive_projections[i] = reduced_descriptor_projections[i]
        else:
            negative_projections[i] = -reduced_descriptor_projections[i]

    final_cartridge_projections = np.concatenate([positive_projections, negative_projections])

    # Step 6: Normalize the final projection vector
    cartridge_percentages = normalize(final_cartridge_projections.reshape(1, -1), norm='l1') * 100

    # Step 7: Calculate approximation accuracy
    variance_explained = np.sum(pca.explained_variance_ratio_) * 100

    return {
        "cartridge_activations": cartridge_percentages[0].tolist(),
        "approximation_accuracy": variance_explained
    }


def get_smiles_from_name(name):
    """Retrieve SMILES string for a compound name."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
        else:
            print(f"Failed to retrieve SMILES for {name}")
            return None
    except Exception as e:
        print(f"Error retrieving SMILES for {name}: {e}")
        return None
