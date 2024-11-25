from flask import Flask, request, jsonify, send_from_directory
from pipeline import query_openai_api, calculate_cartridge_activations
from flask_cors import CORS
import json
from utils import clean_openai_response
import os

print(f"Using query_openai_api from {query_openai_api.__module__}")
app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})



# Predefined SDF file paths for cartridges
SDF_FILES = [
    "smell_UI_CAS_SDF/1.112-37-8.sdf",
    "smell_UI_CAS_SDF/2.105-87-3.sdf",
    "smell_UI_CAS_SDF/3.2277-19-2.sdf",
    "smell_UI_CAS_SDF/4.8007-80-5.sdf",
    "smell_UI_CAS_SDF/5.25634-93-9.sdf",
    "smell_UI_CAS_SDF/6.106-29-6.sdf",
    "smell_UI_CAS_SDF/7.1205-17-0.sdf",
    "smell_UI_CAS_SDF/8.2244-16-8.sdf",
    "smell_UI_CAS_SDF/9.89-43-0.sdf",
    "smell_UI_CAS_SDF/10.31906-04-4.sdf",
]

# Route to serve SDF files
@app.route('/smell_UI_CAS_SDF/<path:filename>', methods=['GET'])
def serve_sdf_file(filename):
    # Define the absolute path to the SDF directory
    sdf_directory = os.path.join(os.getcwd(), 'data/smell_UI_CAS_SDF')
    return send_from_directory(sdf_directory, filename)

@app.route("/run_pipeline", methods=["POST"])
def run_pipeline():
    # Add a log to confirm the request is received
    print("###Request received:", request.json)
    try:
        data = request.json
        print(f"###Incoming request data: {data}")  # Log the incoming request

        smell_description = data.get("smell_description")
        if not smell_description:
            print("###Missing 'smell_description'")  # Log missing key
            return jsonify({"error": "Missing 'smell_description' in request"}), 400

        # Mock query response
        compound_mix_raw = query_openai_api(smell_description)
        print(f"###openAI query response: {compound_mix_raw}")  # Log the response

        # Clean response
        compound_mix = clean_openai_response(compound_mix_raw)
        print(f"###jsonified response: {compound_mix}")

        if not compound_mix or "odorants" not in compound_mix:
            print("###Invalid compound mix format")  # Log invalid mix
            return jsonify({"error": "Invalid compound mix format"}), 500
        
        print(f"###calling cartridge_activations with: {compound_mix['odorants']}")
        # Calculate cartridge activations and approximation accuracy
        cartridge_data = calculate_cartridge_activations(compound_mix["odorants"])
        cartridge_activations = cartridge_data["cartridge_activations"]
        approximation_accuracy = cartridge_data["approximation_accuracy"]

        # print(f"###Cartridge activations returned: {cartridge_activations}")  # Log activations

        # Prepare cartridge data with SDF file paths
        BASE_SDF_URL = "https://laughing-adventure-9qggv74gp54395xp-5000.app.github.dev/smell_UI_CAS_SDF/"
        cartridge_response = [
            {"id": i + 1, "name": f"Cartridge {i + 1}", "activation": activation, "sdfFile": f"{BASE_SDF_URL}{os.path.basename(SDF_FILES[i])}"}
            for i, activation in enumerate(cartridge_activations)
        ]
        return jsonify({
            "compoundMix": compound_mix,
            "cartridge_activations": cartridge_response,
            "approximation_accuracy": approximation_accuracy
        }), 200
    except Exception as e:
        print(f"###Unexpected error: {e}")  # Log unexpected errors
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
