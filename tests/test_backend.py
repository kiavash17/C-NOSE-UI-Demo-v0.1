import pytest
from unittest.mock import patch  # Import mocking tools
from backend.app import app  # Import the Flask app

# Mock the query_openai_api function to isolate the test
@patch("backend.pipeline.query_openai_api", return_value='{"odorants": [{"name": "compound1", "percentage": 50}]}')
def test_run_pipeline(mock_query_openai_api):
    # Create a test client for the Flask app
    client = app.test_client()
    
    # Send a POST request to the endpoint
    response = client.post(
        "/run_pipeline",
        json={"smell_description": "ocean breeze"},  # Test input
        headers={"Content-Type": "application/json"}  # Explicitly set Content-Type
    )
    
    # Assert the status code is 200 (success)
    assert response.status_code == 200

    # Parse the response data
    data = response.get_json()
    
    # Assert that the expected keys exist in the response
    assert "compoundMix" in data
    assert "cartridge_activations" in data
