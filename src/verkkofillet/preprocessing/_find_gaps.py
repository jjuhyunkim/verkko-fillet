# Ensure that 'gaps' column exists
import pandas as pd
import copy 

def find_elements_with_brackets(input_list):
    # Use list comprehension to find indices of elements that start with "[" and end with "]"
    indices = [idx for idx, element in enumerate(input_list) if element.startswith("[") and element.endswith("]")]
    return indices

def findGaps(ob_ori):
    obj = copy.deepcopy(ob_ori)
    
    gapDb = pd.DataFrame()
    if 'gaps' not in obj.paths.columns:
        obj.paths['gaps'] = None  # You can initialize with None or an empty list
    
    # Iterate through all rows in the DataFrame
    for contig in range(0, obj.paths.shape[0]):
       
        # Check if 'path' is not NaN or empty before processing
        input_list = obj.paths.loc[contig, "path"]
        contig_name = obj.paths.loc[contig, "name"]
        if pd.isna(input_list) or not input_list:
            continue  # Skip this row if 'path' is NaN or empty
        
        # Split the path string into a list by commas
        input_list = input_list.split(',')
        
        # Find indices of all elements that start with "[" and end with "]"
        indices = find_elements_with_brackets(input_list)
        # print(f"Indices for row {contig}: {indices}")  # Debugging line
        
        gap_list = []
        
        # Ensure there's no out-of-bounds error by checking the indices
        for idx in indices:
            if 0 < idx < len(input_list) - 1:  # Ensure that idx-1 and idx+1 are valid
                gap_list.append([input_list[idx-1], input_list[idx], input_list[idx+1]])
            else:
                print(f"Out of bounds index at {idx} for row {contig}")  # Debugging line
        
        # Debug print the gap_list
        # print(f"gap_list for row {contig}: {gap_list}")
        
        # Assign the result to the 'gaps' column
        # Use .at for a single row assignment
        gapDb_tmp = pd.DataFrame({'name' : contig_name,
                                  'gaps' : gap_list})
        gapDb = pd.concat([gapDb, gapDb_tmp]) 
    
    # Check the result (optional)
    # print(obj.paths.head())
    gapDb['gapId'] = ["gapid_" + str(i) for i in range(0, gapDb.shape[0])]
    print(str(gapDb.shape[0]) + ' gaps were found -> obj.gaps')
    obj.gaps = gapDb.copy()
    obj.gaps["notes"] = ""
    obj.gaps["fixedPath"] = ""
    obj.gaps["done"] = ""
    obj.gaps["startMatch"] = ""
    obj.gaps["endMatch"] = ""
    obj.gaps["finalGaf"] = ""
    obj.gaps = obj.gaps.loc[:,['gapId','name','gaps','notes','fixedPath','startMatch','endMatch','finalGaf','done']]
    return obj