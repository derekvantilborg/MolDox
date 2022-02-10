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

- [conda](https://anaconda.org/) | Conda must be installed first
- [meeko v0.2](https://github.com/forlilab/Meeko) |```pip install meeko==0.2```
- [Jupyter](https://jupyter.org/) |```pip install jupyter```


## Installation
***


It is advised to create separate conda environment

```
conda create -n moldox python=3.8
conda activate moldox
pip install meeko==0.2
pip install jupyter
```

### Conda installation

MolDox can be simple be installed as

```conda install -c derekvantilborg moldox```

### Manual installation
In special cases, you can install it manually as well.

First clone the git repository:

```git clone https://github.com/derekvantilborg/MolDox```

Install pymol, openbabel, rdkit, MDAnalysis, pdbfixer, vina, py3Dmol, meeko, and jupyter
```
conda install -c conda-forge -c schrodinger pymol-bundle
conda install -c conda-forge openbabel
conda install -c conda-forge rdkit
conda install -c conda-forge MDAnalysis
conda install -c conda-forge pdbfixer
conda install -c conda-forge vina
conda install -c conda-forge py3Dmol
```

## Getting started
***

The easiest way to get started is opening up the included Jupyter Notebook.
Simply type ```jupyter-notebook``` in the terminal and open ```moldox.ipynb```.

An example of this notebook can be found [here](https://github.com/derekvantilborg/MolDox/blob/main/moldox.ipynb). It features fetching
a pdb + ligand from the internet and re-docking it.





