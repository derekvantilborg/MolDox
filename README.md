# MolDox
### Easy to use docking tool based on autodock Vina.
**MolDox takes care of prepping your receptor and ligand(s) from various input formats, defining the 
docking box dimensions, fast docking, and visualizing the results without the use of external executables.**


Inspired by https://github.com/AngelRuizMoreno/Jupyter_Dock

![MolDox logo](img/MolDox.png?raw=true "Title")

***

## Installation
Create conda environment
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





