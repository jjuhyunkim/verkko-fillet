import os
import pandas as pd

def getQV(obj, qvFile="kmer/assembly.qv_cal.qv"):
    """
    Reads a QV (Quality Value) file, parses it, and attaches the resulting DataFrame to the provided object.
    
    Parameters
    ----------
    obj
        The object to which the DataFrame will be attached.
    qvFile
        The path to the QV file. Default is "kmer/assembly.qv_cal.qv".
    
    Returns
    -------
    obj : object
        The object with the QV DataFrame attached
    """
    # Construct the full path to the QV file
    qvFilePath = os.path.abspath(qvFile)

    # Check if the file exists
    if not os.path.exists(qvFilePath):
        raise FileNotFoundError(f"QV file not found: {qvFilePath}")
    
    try:
        # Read the file into a DataFrame
        qvTab = pd.read_csv(qvFilePath, header=None, sep='\t')
        
        # Assign column names
        qvTab.columns = ['asmName', 'nKmer_uniq_asm', 'nKmer_total', 'QV', 'ErrorRate']
        
        # Attach the DataFrame to the object
        obj.qv = qvTab
    
    except Exception as e:
        raise ValueError(f"Error reading or processing the QV file: {e}")
    
    return obj