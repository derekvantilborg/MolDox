
from docking.dock import MolDox
from pymol import cmd
from rdkit import Chem


# cmd.fetch(code='1AZ8',type='pdb1')
# cmd.select(name='Prot',selection='polymer.protein')
# cmd.select(name='Lig',selection='organic')
# cmd.save(filename='test_output/1AZ8_clean.pdb',format='pdb',selection='Prot')
# cmd.save(filename='test_output/1AZ8_lig.mol2',format='mol2',selection='Lig')
# cmd.delete('all')

# cmd.fetch(code='1AZ8',type='pdb1')
cmd.load("example/4ypq.pdb")
cmd.select(name='Prot', selection='polymer.protein')
cmd.select(name='Lig', selection='organic')
cmd.save(filename='test_output2/4ypq_clean.pdb',format='pdb',selection='Prot')
cmd.save(filename='test_output2/4ypq_lig.mol2',format='mol2',selection='Lig')
cmd.delete('all')

dox = MolDox(output_dir='test_output2')
dox.prep_receptor(receptor='test_output2/4ypq_clean.pdb',
                  ligand='test_output2/4ypq_lig.mol2',
                  ph=7.4, keep_water=True, renumber=True)

# dox.prep_reference_ligand('test_output2/4ypq_lig.mol2')

dox.prep_ligands('example/RORytligand_testDerek.sdf', format='sdf')

dox.autodock_all()


dox.docking_results

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
# fix_protein(filename='1AZ8_clean.pdb', addHs_pH=7.4, try_renumberResidues=True, output='1AZ8_clean_H.pdb')
# from docking.prepare_receptor import compute_partial_charges
# compute_partial_charges('1AZ8_clean_H.pdb', 'pdb', pdbqt_output='1AZ8_clean_H_prepped.pdbqt', verbose=False)

prep_receptor('1AZ8_clean.pdb', ph=7.4, keep_water=True, renumber=True, verbose=True)

# Ligand
sanitize_mol2('1AZ8_lig.mol2', '1AZ8_lig_H.mol2')
prepare_ligand('1AZ8_lig_H.mol2', '1AZ8_lig_H.pdbqt')


#################
# Define Box    #
#################

# cmd.load(filename='1AZ8_clean_H.pdb', format='pdb', object='prot')
# cmd.load(filename='1AZ8_lig_H.mol2', format='mol2', object='lig')
#
# center, size = getbox(selection='lig', extending=5.0)
#
# cmd.delete('all')

center, size = find_box(receptor='1AZ8_clean_H.pdb',
                        ligand='1AZ8_lig_H.mol2',
                        receptor_format='pdb',
                        ligand_format='mol2', box_extension=5)


#################
# Autodock Vina #
#################

v = Vina(sf_name='vina')

v.set_receptor('1AZ8_clean_H.pdbqt')

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

import random
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]


#

poses = Chem.SDMolSupplier('test_output/1AZ8_lig_vina_out.sdf',True)
for p in list(poses):
    pose_1 = Chem.MolToMolBlock(p)
    # print(p.GetProp('_Name'), 'Score: {}'.format(p.GetProp('minimizedAffinity')))

    print(p.GetProp('_Name'), p.GetProp('Score'))
    # p.GetProp('')



cmd.fetch(code='1X1R',type='pdb1')
cmd.select(name='Prot',selection='polymer.protein')
cmd.select(name='GDP',selection='organic')
cmd.save(filename='1X1R_clean.pdb',format='pdb',selection='Prot')
cmd.save(filename='1X1R_GDP.mol2',format='mol2',selection='GDP')
cmd.delete('all')




