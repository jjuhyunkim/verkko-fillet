# Installation

Install `verkko-fillet` via PyPI by running:

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
