from flask import Flask, request, jsonify
from pipeline import query_openai_api, calculate_cartridge_activations
from flask_cors import CORS
import json

app = Flask(__name__)
CORS(app)  # Enable Cross-Origin Resource Sharing

@app.route("/run_pipeline", methods=["POST"])
def run_pipeline():
    try:
        # Step 1: Get the incoming data (smell description)
        data = request.json
        print(f"Incoming request data: {data}")
        
        smell_description = data.get("smell_description")
        if not smell_description:
            return jsonify({"error": "Missing 'smell_description' in request"}), 400

        # Step 2: Query OpenAI for compound mix
        compound_mix_raw = query_openai_api(smell_description)
        if "error" in compound_mix_raw:
            return jsonify(compound_mix_raw), 400

        # Step 3: Calculate cartridge activations (handle response and process)
        compound_mix = json.loads(compound_mix_raw)
        if not compound_mix or "odorants" not in compound_mix:
            return jsonify({"error": "Invalid compound mix format"}), 500
        
        # Step 4: Calculate cartridge activations based on the odorants returned
        cartridge_activations = calculate_cartridge_activations(compound_mix["odorants"])

        # Step 5: Return response
        return jsonify({"compoundMix": compound_mix, "cartridge_activations": cartridge_activations}), 200

    except Exception as e:
        print(f"Unexpected error: {e}")
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
