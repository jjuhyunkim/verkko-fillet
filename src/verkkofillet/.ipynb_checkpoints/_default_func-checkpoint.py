import os 
import pandas as pd
from datetime import datetime


def addHistory(obj, activity, function):
    # Initialize an empty DataFrame for history
    history = obj.history.copy()
    # Get the current timestamp
    now = datetime.now()
    
    # Create a new history entry
    new_history = pd.DataFrame({'timestamp': [now], 'activity': [activity], 'function' : [function]})
    
    # Check if 'history' is empty before concatenating
    if not history.empty:
        history = pd.concat([history, new_history], ignore_index=True)
    else:
        # If the history DataFrame is empty, directly assign the new history
        history = new_history
    
    obj.history = history
    # Display the updated history DataFrame
    return obj


def check_user_input(user_input, default_value):
    """
    Checks if the user input matches the default value. 
    Returns the default value if matched, otherwise returns the user input.
    """
    return user_input if user_input != default_value else default_value


def print_directory_tree(base_path, max_depth=1, prefix="", is_root=True):
    """
    Prints a directory tree structure from and including the base directory.
    
    :param base_path: Path to the directory to scan.
    :param max_depth: Maximum depth to explore.
    :param prefix: Prefix for tree structure (used internally for recursion).
    :param is_root: Flag to indicate if this is the root directory.
    """
    if is_root:
        print(base_path)  # Print the base directory at the top
    
    if max_depth < 1:
        return

    try:
        # List directory contents
        entries = os.listdir(base_path)
    except PermissionError:
        print(f"{prefix}[Permission Denied]")
        return

    for index, entry in enumerate(entries):
        entry_path = os.path.join(base_path, entry)
        is_last = (index == len(entries) - 1)

        # Use "├── " or "└── " to show branches in the tree
        connector = "└── " if is_last else "├── "
        print(f"{prefix}{connector}{entry}")

        # Recurse into subdirectories
        if os.path.isdir(entry_path):
            # Add prefix for subdirectories
            sub_prefix = "    " if is_last else "│   "
            print_directory_tree(entry_path, max_depth - 1, prefix + sub_prefix, is_root=False)