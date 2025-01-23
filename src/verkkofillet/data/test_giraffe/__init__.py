import os

# Function to get the path of a script
def get_script_path(script_name):
    return os.path.join(os.path.dirname(__file__), script_name)