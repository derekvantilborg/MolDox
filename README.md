# MolDox
***
### Easy to use docking tool based on autodock Vina.
**MolDox takes care of prepping your receptor and ligand(s) from various input formats, defining the 
docking box dimensions, fast docking, and visualizing the results without the use of external executables.**


Inspired by https://github.com/AngelRuizMoreno/Jupyter_Dock

![MolDox logo](img/MolDox.png?raw=true "Title")


## Features
***

- Easy and fast docking of ligands using Vina Autodock 
- Automatic docking box dimensions detection
- Automatic receptor and ligand prepping
- Interactive visualization of docking results
- Simulation of Protein-ligand interactions 

## Requirements
***
MolDox currently supports Python 3.8

- [vina](https://vina.scripps.edu/)
- [pymol](https://anaconda.org/schrodinger/pymol-bundle)
- [openbabel](https://anaconda.org/openbabel/openbabel)
- [pdbfixer](https://anaconda.org/omnia/pdbfixer)
- [MDAnalysis](https://www.mdanalysis.org/pages/installation_quick_start/)
- [rdkit](https://www.rdkit.org/)
- [py3Dmol](https://github.com/avirshup/py3dmol)
- [meeko](https://github.com/forlilab/Meeko)
- [ProLIF](https://github.com/chemosim-lab/ProLIF)

## Installation
***
MolDox can *NOT YET* be installed as

```pip install git+https://github.com/derekvantilborg/MolDox.git```

### Manual installation
In the mean time you can install it manually:

```git clone https://github.com/derekvantilborg/MolDox```

Create a conda environment

```
conda create -n docking python=3.8
conda activate moldox
```
Install pymol, openbabel, rdkit, MDAnalysis, pdbfixer
```
conda install -c conda-forge -c schrodinger pymol-bundle
conda install -c conda-forge openbabel
conda install -c conda-forge rdkit
conda install -c conda-forge MDAnalysis
conda install -c conda-forge MDAnalysisTests
conda install -c conda-forge pdbfixer
```
Pip install vina, py3Dmol, meeko, ProLIF, and jupyter
```
pip install vina
pip install py3Dmol
pip install meeko
pip install cyton  # required for ProLIF
pip install git+https://github.com/chemosim-lab/ProLIF.git
pip install jupyter
```

## Getting started
***

A tutorial can be found [here](https://github.com/derekvantilborg/MolDox/blob/main/moldox.ipynb) that features fetching
a pdb + ligand from the internet and re-docking it.





