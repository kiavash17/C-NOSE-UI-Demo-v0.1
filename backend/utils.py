import requests
import json

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
    
def clean_openai_response(response: str) -> dict:
    """
    Cleans and parses the response from OpenAI, handling common formatting issues
    such as extra backticks and 'json' prefixes.

    Args:
        response (str): Raw response string from OpenAI.

    Returns:
        dict: Parsed JSON response as a dictionary.

    Raises:
        ValueError: If the response cannot be cleaned or parsed as valid JSON.
    """
    try:
        # Strip backticks and whitespace
        cleaned_response = response.strip("`").strip()

        # Remove "json\n" prefix if it exists
        if cleaned_response.startswith("json"):
            cleaned_response = cleaned_response[5:].strip()

        # Parse the cleaned string as JSON
        return json.loads(cleaned_response)

    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON response: {e}")