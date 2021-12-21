from meeko import MoleculePreparation
from meeko import obutils
from openbabel import pybel
from openbabel import openbabel
import os

from rdkit import Chem
from rdkit.Chem import Draw
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import py3Dmol
from pymol import cmd
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
import random, math
import numpy as np
from vina import Vina

from MDAnalysis.coordinates import PDB

import prolif as plf
from prolif.plotting.network import LigNetwork

from docking.prepare import fix_protein, sanitize_mol2, getbox, prepare_ligand, hydrogenate_and_compute_partial_charges, pdbqt_to_sdf

#################
# Fetch pdb     #
#################

cmd.fetch(code='1AZ8',type='pdb1')
cmd.select(name='Prot',selection='polymer.protein')
cmd.select(name='Lig',selection='organic')
cmd.save(filename='1AZ8_clean.pdb',format='pdb',selection='Prot')
cmd.save(filename='1AZ8_lig.mol2',format='mol2',selection='Lig')
cmd.delete('all')


#################
# Sanitization  #
#################

# Protein
fix_protein(filename='1AZ8_clean.pdb', addHs_pH=7.4, try_renumberResidues=True, output='1AZ8_clean_H.pdb')
hydrogenate_and_compute_partial_charges('1AZ8_clean_H.pdb', 'pdb', hyd_output=None,
                                        pdbqt_output='1AZ8_clean_H_prepped.pdbqt',
                                        protein=True, verbose=False)

# Ligand
sanitize_mol2('1AZ8_lig.mol2', '1AZ8_lig_H.mol2')
prepare_ligand('1AZ8_lig_H.mol2', '1AZ8_lig_H.pdbqt')


#################
# Define Box    #
#################

cmd.load(filename='1AZ8_clean_H.pdb', format='pdb', object='prot')
cmd.load(filename='1AZ8_lig_H.mol2', format='mol2', object='lig')

center, size = getbox(selection='lig', extending=5.0)

cmd.delete('all')


#################
# Autodock Vina #
#################

v = Vina(sf_name='vina')

v.set_receptor('1AZ8_clean_H_prepped.pdbqt')

v.set_ligand_from_file('1AZ8_lig_H.pdbqt')

v.compute_vina_maps(center=[center['center_x'], center['center_y'], center['center_z']],
                    box_size=[size['size_x'], size['size_y'], size['size_z']])

'''
# Score the current pose
energy = v.score()
print('Score before minimization: %.3f (kcal/mol)' % energy[0])

# Minimized locally the current pose
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('1iep_ligand_minimized.pdbqt', overwrite=True)
'''

# Dock the ligand
v.dock(exhaustiveness=10, n_poses=10)
v.write_poses('1AZ8_lig_vina_out.pdbqt', n_poses=10, overwrite=True)

###########################
# Docking results to SDF #
###########################

pdbqt_to_sdf(pdbqt_file='1AZ8_lig_vina_out.pdbqt', output='1AZ8_lig_vina_out.sdf')

####

results = Chem.SDMolSupplier('1AZ8_lig_vina_out.sdf')
p = Chem.MolToMolBlock(results[0], False)



#####

fix_protein(filename='1AZ8_clean.pdb', addHs_pH=7.4, try_renumberResidues=True, output='1AZ8_clean_H_fix.pdb')

# load protein
prot = mda.Universe("1AZ8_clean_H_fix.pdb")
prot = plf.Molecule.from_mda(prot)
prot.n_residues

# load ligands
lig_suppl = list(plf.sdf_supplier('1AZ8_lig_vina_out.sdf'))
# generate fingerprint
fp = plf.Fingerprint()
fp.run_from_iterable(lig_suppl, prot)
results_df = fp.to_dataframe(return_atoms=True)
results_df


net = LigNetwork.from_ifp(results_df, lig_suppl[0], kind="frame", frame=0, rotation=270)
net.display()


