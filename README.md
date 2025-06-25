# verkko-fillet : post verkko graph and assembly cleaning in Python

verkko-fillet is an easy-to-use toolkit for cleaning [Verkko](https://github.com/marbl/verkko) assemblies. It includes tools for generating the Verkko-Fillet object, performing assembly quality checks, identifying gaps, assigning chromosomes, searching for ONT reads to help resolve gaps, filling gaps, and generating the final Rukki path (in a GAF-like format) for future Verkko CNS runs.

This Python-based implementation streamlines the entire process, starting right after the Verkko assembly is completed and preparing for the CNS run.


### Installation

dependencies : 
* [gfacpp](https://github.com/snurk/gfacpp)

Using `pip` is recommended. [link](https://pypi.org/project/verkkofillet/)


The default name of the Mamba or Conda environment is `verkko-fillet`. If you want to use a different name, please update the name field in the `environment.yaml` file before proceeding.

```
# Generate a mamba or conda environment and install dependencies.
mamba create -n verkko-fillet -f vf_environment.Jun252025.yml -vvv --dry-run --channel-priority flexible

# Once you’ve checked the generated environment, please re-run this without the --dry-run option. It may take some time.
mamba create -n verkko-fillet -f vf_environment.Jun252025.yml -vvv --channel-priority flexible

# Acivate the environment.
mamba activate verkko-fillet # or the name you desired

# Add python jupyter kernel.
python -m ipykernel install --user --name verkko-fillet --display-name verkko-fillet
pip install verkkofillet
```