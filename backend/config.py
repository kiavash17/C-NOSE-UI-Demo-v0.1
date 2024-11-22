from dotenv import load_dotenv
import os

# Load environment variables from .env
load_dotenv()

# Access the OpenAI API key
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")