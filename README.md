# verkko-fillet : post verkko graph and assembly cleaning in Python

verkko-fillet is an easy-to-use toolkit for cleaning [Verkko](https://github.com/marbl/verkko) assemblies. It includes tools for generating the Verkko-Fillet object, performing assembly quality checks, identifying gaps, assigning chromosomes, searching for ONT reads to help resolve gaps, filling gaps, and generating the final Rukki path (in a GAF-like format) for future Verkko CNS runs.

This Python-based implementation streamlines the entire process, starting right after the Verkko assembly is completed and preparing for the CNS run.

dependencies : 
* gfacpp

### Installation
Using `pip` is recommended. [link](https://pypi.org/project/verkkofillet/)
```
mamba create -n verkko-fillet python=3.9
mamba activate verkko-fillet
mamba env update --file environment.yaml
pip install verkkofillet==0.1.1
```
