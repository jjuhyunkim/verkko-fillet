import pandas as pd
import logging
from .._default_func import addHistory
# Configure logging
logging.basicConfig(level=logging.INFO)

def progress_bar(current, total):
    """
    Displays a progress bar in the console.
    Args:
        current (int): Current progress.
        total (int): Total progress.
    """
    progress = int((current / total) * 50)
    bar = "[" + "=" * progress + " " * (50 - progress) + "]"
    print(f"\r{bar} {current}/{total} gaps filled", end="")
    print(" ")

def checkGapFilling(obj):
    """
    This function checks and prints the number of filled gaps in the 'gap' DataFrame
    and shows the progress bar for gap filling.

    Args:
        obj: An object that contains the 'gap' DataFrame in obj.gaps.
    """
    total = obj.gaps.shape[0]  # Total number of gaps
    gap = obj.gaps  # Assuming gap is the DataFrame containing gap information

    gap['finalGaf'] = gap['finalGaf'].str.replace('<', '&lt;').str.replace('>', '&gt;')
    gap['done'] = gap['finalGaf'].apply(lambda x: "✅" if x else "")
    # Count the number of non-empty 'finalGaf' entries
    current = gap['finalGaf'].apply(lambda x: pd.notna(x) and x != "").sum()
    
    # Print the current and total number of filled gaps
    # print(f"Number of filled gaps: {current} of total gaps: {total}")

    # Call the progress_bar function to show the filling progress
    progress_bar(current, total)
    
    return gap

def transform_path(elements):
    """
    Transforms elements of the path for gap filling.
    Args:
        elements (list): List of path elements.
    Returns:
        list: Transformed path elements.
    """
    return [
        (">" + elem[:-1] if elem.endswith("+") else "<" + elem[:-1]) if not elem.startswith("[") else elem
        for elem in elements
    ]

def check_match(gap_value, element, position):
    """
    Checks if a specific gap matches the given element.
    Args:
        gap_value (str): The gap value from the DataFrame.
        element (str): The element to match.
        position (int): The position in the gap (0 for start, 2 for end).
    Returns:
        str: "match" if matches, else "notMatch".
    """
    return "match" if gap_value[position] == element else "notMatch"

def fillGaps(obj, gapId, final_path):
    """
    Fills gaps for a specific gapId, updates the 'fixedPath', 'startMatch', 'endMatch', and 'finalGaf' columns.
    
    Args:
        obj: The object containing gap data.
        gapId: The specific gap ID to fill.
        final_path: The path for gap filling. If empty, the relevant columns are filled with 'NA'.
    """
    gap = obj.gaps  # The DataFrame containing gap data

    # Ensure the gapId exists
    if gapId not in gap['gapId'].values:
        raise ValueError(f"gapId {gapId} not found in the DataFrame.")

    # Handle empty final_path
    if final_path == "":
        gap.loc[gap['gapId'] == gapId, ['fixedPath', 'startMatch', 'endMatch', 'finalGaf']] = ""
        print(f"gapId {gapId}: 'final_path' is empty. Other columns have been reset to 'NA'.")
    else:
        # Update the 'fixedPath' column for the matching gapId
        gap.loc[gap['gapId'] == gapId, 'fixedPath'] = final_path

        elements = final_path.replace(" ", "").split(",")
        modified_elements = transform_path(elements)
        modified_path = "".join(modified_elements)
        print(f"final path : {modified_path}")

        # Update the 'finalGaf' column for the matching gapId
        gap.loc[gap['gapId'] == gapId, 'finalGaf'] = modified_path

        # Retrieve the matching row for further updates
        gap_row = gap.loc[gap['gapId'] == gapId].iloc[0]

        # Check the direction and update 'startMatch' and 'endMatch'
        gap.loc[gap['gapId'] == gapId, 'startMatch'] = check_match(gap_row.gaps, elements[0], 0)
        gap.loc[gap['gapId'] == gapId, 'endMatch'] = check_match(gap_row.gaps, elements[-1], 2)

        print(f"Updated gapId {gapId}!")
        print(" ")
        if check_match(gap_row.gaps, elements[0], 0) == "match" :
            print("✅ The start node and its direction match the original node.")
        else :
            print("❌ The start node and its direction do not match the original node.")
        
        if check_match(gap_row.gaps, elements[-1], 2) == "match" :
            print("✅ The start node and its direction match the original node.")
        else :
            print("❌ The start node and its direction do not match the original node.")
        
    # Count remaining empty strings or 'NA' in 'finalGaf
    obj.gaps = gap
    
    obj = addHistory(obj, f"{gapId} filled with {final_path}", 'fillGaps')
    # Show progress after each gap filled
    checkGapFilling(obj)
    
    # Return the updated object
    
    return obj




# Reset the index of the 'gap' DataFrame
def writeFixedGaf(obj, rukki = "8-hicPipeline/rukki.paths.gaf", save = "final_rukki_fixed.paths.gaf"):
    print("Reading original rukki path from "+rukki)
    ori_rukki = pd.read_csv(rukki, sep ='\t', header = 0)
    
    gap = obj.gaps.reset_index()
    
    # Iterate through each row of the 'gap' DataFrame
    for num in range(0, gap.shape[0]):  # Use shape[0] for rows, not shape[1]
        finalGaf = gap.loc[num, 'finalGaf']
        
        # If finalGaf is an empty string, stop the iteration
        if finalGaf == "":
            # print("Stopping as finalGaf is empty.")
            break  # Exit the loop if finalGaf is empty
        
        # Extract relevant values
        contig = gap.loc[num, 'name']
        ori_gap = gap.loc[num, 'gaps']
        
        # Fetch the 'path' for the given contig from 'ori_rukki'
        ori_path_series = ori_rukki.loc[ori_rukki['name'] == contig, 'path']
        
        # Check if ori_path_series is not empty to avoid IndexError
        if ori_path_series.empty:
            print(f"Warning: No path found for contig '{contig}'")
            continue  # Skip this iteration if no path is found
    
        ori_path = ori_path_series.tolist()[0]
    
        # Transform the ori_gap using the transform_path function
        modified_elements = transform_path(ori_gap)
        modified_path = "".join(modified_elements)
        
        # Replace the modified path with the finalGaf, handling HTML escape sequences
        fixedGaf = ori_path.replace(modified_path, finalGaf).replace('&lt;', '<').replace('&gt;', '>')
    
        # Update the ori_rukki DataFrame with the new fixedGaf value
        ori_rukki.loc[ori_rukki['name'] == contig, 'path'] = fixedGaf
    
        print(f"Updated path for contig '{contig}': {fixedGaf}")
    
    # After the loop, ori_rukki will be updated, and you can do something with it, like printing
    # print(ori_rukki)
    ori_rukki.to_csv(save, index=False, sep = '\t')
    obj = addHistory(obj, f"Final rukki path was saved as {save}", 'writeFixedGaf')
    print("Writing fixed rukki path to "+save)