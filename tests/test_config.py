import os
from dotenv import load_dotenv

def test_openai_api_key():
    # Ensure the .env file is loaded
    load_dotenv()
    api_key = os.getenv("OPENAI_API_KEY")

    # Assert the key is set
    assert api_key is not None, "API key is not set!"
    assert len(api_key) > 0, "API key is empty!"
