import pandas as pd
import logging
import re
import copy
from .._default_func import addHistory
# Configure logging
logging.basicConfig(level=logging.INFO)

def path_to_gaf(input_string):
    # Split the input string by commas
    items = input_string.split(',')
    
    # Process each item
    result = []
    for item in items:
        if item.endswith('+'):
            result.append('>' + item[:-1])
        elif item.endswith('-'):
            result.append('<' + item[:-1])
        elif item.endswith(']'):
            result.append(item)
        else :
            error = "Invalid input string: " + input_string
            print(error)
            return error

    
    # Join the processed items into a single string
    return ''.join(result)


def progress_bar(current, total):
    """
    Displays a progress bar in the console.
    
    Parameters
    ----------
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

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    
    Returns
    ----------
        The updated 'gap' DataFrame.
    """
    obj_sub = copy.deepcopy(obj)
    total = obj_sub.gaps.shape[0]  # Total number of gaps
    gap = obj_sub.gaps  # Assuming gap is the DataFrame containing gap information

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

    Parameters
    ----------
    elements
        A list of elements in the path.
    
    Returns
    -------
    list
        A list of transformed elements.
    """
    return [
        (">" + elem[:-1] if elem.endswith("+") else "<" + elem[:-1]) if not elem.startswith("[") else elem
        for elem in elements
    ]

def check_match(gap_value, element, position):
    """
    Checks if a specific gap matches the given element.
    
    Parameters
    ----------
    gap_value
        The gap value from the DataFrame.
    element
        The element to match.
    position
        The position in the gap (0 for start, 2 for end).
    
    Returns
    -------
        "match" if matches, else "notMatch".
    """
    return "match" if gap_value[position] == element else "notMatch"

def fillGaps(obj, *, gapId, final_path):
    """
    Fills gaps for a specific gapId, updates the 'fixedPath', 'startMatch', 'endMatch', and 'finalGaf' columns.
    
    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    gapId
        The identifier for the gap.
    final_path
        The final path to fill the gap.

    Returns
    -------
    obj : object
        The updated verkko fillet object.
    """
    obj_sub = copy.deepcopy(obj)
    gap = obj_sub.gaps  # The DataFrame containing gap data

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
    obj_sub.gaps = gap
    
    obj_sub = addHistory(obj_sub, f"{gapId} filled with {final_path}", 'fillGaps')
    # Show progress after each gap filled
    checkGapFilling(obj_sub)
    
    # Return the updated object
    
    return obj_sub

# Reset the index of the 'gap' DataFrame


def preprocess_path(path_str):
    path_str = path_str.replace("\[", ",\[")
    path_str = path_str.replace("\]", "\],")
    split_path = re.split(r',', path_str)
    return [p for p in split_path if p.strip()]  # Remove empty elements

def connectContigs(obj, contig, contig_to,  at = "left", gap = "[N5000N:connectContig]", flip = False, fix_path = True):
    """
    Connects two contigs by adding a gap between them.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used. 
    contig
        The name of the contig to connect. 
    contig_to
        The name of the contig to connect to. 
    at
        The position to connect the contig. Default is "left". Other option is "right".
    gap
        The gap to add between the contigs. Default is "[N5000N:connectContig]". 
    flip
        Whether to flip the path from contig. Default is False. If True, the path will be flipped.

    Returns
    -------
        The updated VerkkoFillet object with new gap added.
    """
    obj_sub = copy.deepcopy(obj)
    pathdb = obj_sub.paths.copy()
    gapdb = obj_sub.gaps.copy()

    path1_raw = pathdb.loc[pathdb['name'] == contig_to]["path"].values
    path2_raw = pathdb.loc[pathdb['name'] == contig]["path"].values

    path1_path = preprocess_path(path1_raw[0])
    path2_path = preprocess_path(path2_raw[0])

    gapid_add= f"gapid_{str(gapdb['gapId'].str.replace('gapid_','').astype(int).max()+1)}"

    if flip:
        path2_path = path2_path[::-1]
        path2_path = [s.translate(str.maketrans("+-", "-+")) for s in path2_path]

    if at == "left":
        marker = ["startMarker"]
        fixedPath = path2_path + [gap] + marker
    if at == "right":
        marker = ['endMarker']
        fixedPath = marker + [gap]+  path2_path
        
    fixedPath = ','.join(fixedPath)
    gap_new_line = pd.DataFrame({"gapId": gapid_add, 
                  "name" : [contig_to],
                  "gaps" : marker,
                  "notes" : f"connected {contig} to {contig_to} at {at} with flip {flip}",
                  "fixedPath": fixedPath,
                  "startMatch" : "",
                  "endMatch" : "",
                  "finalGaf" : "",
                  "done" : True})

    gapdb = pd.concat([gapdb,gap_new_line], ignore_index=True)
    obj_sub = addHistory(obj_sub, f"{gapid_add} was created", 'connectContig')
    
    newPath = fixedPath
    pathdb = obj_sub.paths.copy()
    pathdb = pathdb.loc[pathdb["name"] != contig]
    ori_path = pathdb.loc[pathdb["name"] == contig_to, 'path'].values[0]
    
    marker = marker[0]
    if marker == "startMarker":
        pathdb.loc[pathdb["name"] == contig_to, 'path'] = newPath + ori_path
    elif marker == "endMarker":
        pathdb.loc[pathdb["name"] == contig_to, 'path'] = ori_path + newPath
    else:
        print("marker doesn't match")

    gapdb['name'] = gapdb['name'].replace(contig, contig_to)
    
    if fix_path:
        obj_sub.paths = pathdb
    
    obj_sub.gaps = gapdb

    print(f"Connected {contig} to {contig_to} at {at} with flip {flip}")
    print(f"{contig} was merged to {contig_to} in obj.paths")
    print(f"{contig} was replaced with {contig_to} in obj.gaps")
    print(f"New gap was added to obj.gaps with gapId {gapid_add}")
    return obj_sub



def deleteGap(obj, gapId):

    """
    Deletes a gap from the 'gap' DataFrame for a specific gapId.

    Parameters
    ----------
    obj
        An verkko fillet object that contains the 'gap' DataFrame in obj.gaps.
    gapId
        The identifier for the gap to delete.

    Returns
    -------
        The updated verkko fillet object.
    """
    obj_sub = copy.deepcopy(obj)
    if gapId not in obj.gaps['gapId'].values:
        raise ValueError(f"gapId {gapId} not found in the DataFrame.")
    gaps = obj_sub.gaps.copy()
    gaps = gaps.loc[gaps['gapId'] != gapId]
    obj_sub.gaps = gaps
    obj_sub = addHistory(obj_sub, f"{gapId} was removed", 'deleteGap')
    return obj_sub


def writeFixedPaths(obj, save_path = "assembly.fixed.paths.tsv", save_gaf = "assembly.fixed.paths.gaf"):
    """
    Writes the fixed paths to a file.

    Parameters
    ----------
    obj
        The VerkkoFillet object to be used. 
    save_path
        The file path to save the fixed paths. Default is "assembly.fixed.paths.tsv".
    save_gaf
        The file path to save the fixed gaf. Default is "assembly.fixed.paths.gaf".

    Returns
    -------
        The fixed paths and gaf saved to the specified file paths.
    """
    obj_sub = copy.deepcopy(obj)
    gapdb = obj_sub.gaps.copy()
    gapdb = gapdb.reset_index()
    pathdb = obj_sub.paths.copy()

    print("Checking for startMarker and endMarker...")
    print(" ")
    for num in range(len(gapdb)):
        marker = gapdb.loc[num, "gaps"][0]

        if isinstance(marker, str) and re.search(r"startMarker|endMarker", marker):
            print('yes')
            note = gapdb.loc[num, "notes"].split(" ")
            rmnode = note[1]
            remainnode = note[3]

            newPath = gapdb.loc[num, "fixedPath"].replace("startMarker", "").replace("endMarker", "")
            print(newPath)
            # update pathdb
            pathdb = pathdb.loc[pathdb["name"] != rmnode]
            ori_path = pathdb.loc[pathdb["name"] == remainnode, 'path'].values[0]
            print(ori_path)
            # add path to main contig
            if marker == "startMarker":
                pathdb.loc[pathdb["name"] == remainnode, 'path'] = newPath + ori_path
            elif marker == "endMarker":
                pathdb.loc[pathdb["name"] == remainnode, 'path'] = ori_path + newPath
            else:
                print("marker doesn't match")

            gapdb['name'] = gapdb['name'].replace(rmnode, remainnode)

    print("Fixing paths using gap infomation...")
    print(" ")
    for num in range(len(gapdb)):
        marker = gapdb.loc[num, "gaps"][0]

        if isinstance(marker, str) and re.search(r"startMarker|endMarker", marker):
            continue
        if gapdb.loc[num, 'finalGaf'] == "":
            continue
        else:
            fixed_path = gapdb.loc[num, 'fixedPath'].replace(" ", "")
            gap_path = ','.join(gapdb.loc[num, 'gaps'])
            contig = gapdb.loc[num, 'name']

            ori_path = pathdb.loc[pathdb['name'] == contig, 'path'].values[0]
            ori_path = ori_path.replace(gap_path, fixed_path)

            # update path
            pathdb.loc[pathdb['name'] == contig, 'path'] = ori_path

    pathdb = pathdb.reset_index(drop=True)
    pathdb = pathdb.drop(columns = "gaps")
       
    gaf = pathdb.copy()
    gaf['path'] = pathdb['path'].apply(path_to_gaf)

    pathdb.to_csv(save_path, sep = "\t", index = False)
    print(f"Fixed paths were saved to {save_path}")
    
    gaf.to_csv(save_gaf, sep = "\t", index = False)
    print(f"Fixed gaf were saved to {save_gaf}")


