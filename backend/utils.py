import requests

def get_smiles_from_name(name):
    """Fetch SMILES string for a given compound name."""
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
