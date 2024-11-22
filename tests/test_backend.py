import sys
import os
from unittest.mock import patch
import pytest

# Add backend directory to Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../backend')))
from app import app

# Mock the query_openai_api function
@patch("app.query_openai_api", return_value='{"odorants": [{"name": "Vanillin", "percentage": 50}]}')
def test_run_pipeline(mock_query_openai_api):
    # Log to confirm the mock is applied
    print(f"Mock applied: {mock_query_openai_api.called}")  # Should initially be False

    # Create a test client for the Flask app
    client = app.test_client()

    # Send a POST request to the endpoint
    response = client.post(
        "/run_pipeline",
        json={"smell_description": "ocean breeze"},
        headers={"Content-Type": "application/json"}
    )

    # Check if the mock was actually called
    print(f"Mock called: {mock_query_openai_api.called}")  # Should now be True

    # Assert the response status code is 200
    assert response.status_code == 200

    # Parse the response data
    data = response.get_json()

    # Assert the response contains the expected keys
    assert "compoundMix" in data
    assert "cartridge_activations" in data
