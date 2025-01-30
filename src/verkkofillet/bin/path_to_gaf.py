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


# In[2]:


#now = datetime.datetime.now()
#print('TimeStamp : ' + now.strftime('%Y-%m-%d %H:%M:%S'))


# In[6]:


# In[7]:


gaf_file=os.path.realpath(sys.argv[1])
path_file=os.path.realpath(sys.argv[2])

print("ouput GAF file : " + gaf_file)
print("input path.tsv file : " + path_file)


if os.path.exists(gaf_file):
    print(f"GAF file '{gaf_file}' already exists. Exiting the script.")
    sys.exit(1)


# In[8]:


path = pd.read_csv(path_file, header=0, sep='\t')
path.head()


# In[11]:


def transform_string(input_string):
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

# Example usage
#input_string = "utig4-2404-,utig4-2405+,utig4-2687-"
#output_string = transform_string(input_string)
#print(output_string)  # Output: <utig4-2404>utig4-2405<utig4-2687


# In[12]:


gaf = path.copy()


# In[13]:


gaf['path'] = path['path'].apply(transform_string)


# In[14]:


gaf.head()


# In[15]:


print("Writing...")
gaf.to_csv(gaf_file, sep='\t',header=True, index=None)


# In[ ]:

