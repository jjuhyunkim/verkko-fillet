import sys
import shlex
import subprocess
import os

script_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../bin/'))

def findRevChr(obj, size = "large", mashmap_out = "assembly.mashmap.out"):
    print("hi")