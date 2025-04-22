#!/usr/bin/env python
# coding: utf-8

# # Gaf to Path.tsv

# Author : Juhyun Kim 
# Affiliation : NHGRI/NIH

# In[1]:


print("Loading libraries")
#import session_info
import pandas as pd
import csv
import sys
#import datetime
import os
import re


# In[2]:


#now = datetime.datetime.now()
#print('TimeStamp : ' + now.strftime('%Y-%m-%d %H:%M:%S'))


# In[6]:


# In[7]:


gaf_file=os.path.realpath(sys.argv[1])
path_file=os.path.realpath(sys.argv[2])

print("input GAF file : " + gaf_file)
print("output path.tsv file : " + path_file)


if os.path.exists(path_file):
    print(f"path file '{path_file}' already exists. Exiting the script.")
    sys.exit(1)


# In[8]:


gaf = pd.read_csv(gaf_file, header=0, sep='\t')
gaf.head()


# In[11]:


def transform_string(input_string):
    # Split the input string by commas
    # items = input_string.split(',')
    items =   re.findall(r'([><\[][^\]><[]*\]?)', input_string)
    
    
    # Process each item
    result = []
    for item in items:
        if item.startswith('>'):
            result.append(item.lstrip('<>') + '+')
        elif item.startswith('<'):
            result.append(item.lstrip('<>') + '-')
        elif item.endswith(']'):
            result.append(item)
        else :
            error = "Invalid input string: " + item
            print(error)
            return

    
    # Join the processed items into a single string
    return ','.join(result)

# Example usage
#input_string = "utig4-2404-,utig4-2405+,utig4-2687-"
#output_string = transform_string(input_string)
#print(output_string)  # Output: <utig4-2404>utig4-2405<utig4-2687


# In[12]:


path = gaf.copy()


# In[13]:


path['path'] = path['path'].apply(transform_string)


# In[14]:


path.head()


# In[15]:


print("Writing...")
path.to_csv(path_file, sep='\t',header=True, index=None)


# In[ ]:

