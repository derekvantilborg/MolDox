
from meeko import MoleculePreparation
from meeko import obutils
from openbabel import pybel, openbabel
import os
import time

# def sanitize_mol(infile, outfile, overwrite=True):
#     """
#
#     :param infile: (str) pathname to input mol2 file
#     :param outfile: (str) pathname to output mol2 file
#     :param overwrite: (bool) if True, overwrite existing outfile
#     """
#     # read molecule from file
#     mol = [m for m in pybel.readfile(filename=infile, format=infile.split('.')[-1])][0]
#     # Add hydrogens
#     mol.addh()
#     # Write to output file
#     out = pybel.Outputfile(filename=outfile, format=outfile.split('.')[-1], overwrite=overwrite)
#     out.write(mol)
#     out.close()


def prep_mol_from_file(infile, output_pdbqt, hydrate=False, keep_nonpolar_hydrogens=False, pH_value=None):
    mol = obutils.load_molecule_from_file(infile)

    preparator = MoleculePreparation(hydrate=hydrate,
                                     keep_nonpolar_hydrogens=keep_nonpolar_hydrogens,
                                     pH_value=pH_value)
    preparator.prepare(mol)

    preparator.write_pdbqt_file(output_pdbqt)


def hydrogenate_mol_from_file(infile, outfile, overwrite=True):
    """

    :param infile: (str) pathname to input mol2 file
    :param outfile: (str) pathname to output mol2 file
    :param overwrite: (bool) if True, overwrite existing outfile
    """
    # read molecule from file
    mol = [m for m in pybel.readfile(filename=infile, format=infile.split('.')[-1])][0]
    # Add hydrogens
    mol.addh()
    # Write to output file
    out = pybel.Outputfile(filename=outfile, format=outfile.split('.')[-1], overwrite=overwrite)
    out.write(mol)
    out.close()


def prep_ligands(infile, format='sdf', output_dir='.', hydrogenate=True, hydrate=False, keep_nonpolar_hydrogens=False,
                 pH_value=None, make3d=True, remake3d=False, forcefield='mmff94', steps_3d=50, steps_3d_localopt=500):

    # Read file
    mols = [m for m in pybel.readfile(filename=infile, format=format)]

    # Add hydrogens
    if hydrogenate:
        for m in mols:
            m.addh()

    # Make 3D if needed
    for m in mols:
        z_coords = [a.coords[-1] for a in m.atoms]
        if z_coords.count(0) == len(m.atoms):
            print(f'All 0 Z-coordinates of molecule {m}'.strip())
            if make3d:
                print(f'\t3D coordinates are generated with the {forcefield} forcefield ({steps_3d} steps)')
                m.make3D('mmff94', steps_3d)
                if steps_3d_localopt > 0:
                    print(f'\tPerforming local optimization ({steps_3d_localopt} steps)')
                    m.localopt('mmff94', steps_3d_localopt)
        elif remake3d:
            if z_coords.count(0) == len(m.atoms):
                print(f'Existing Z-coordinates, but new 3D coordinates are computed anyways for {m}'.strip())
            print(f'\t3D coordinates are generated with the {forcefield} forcefield ({steps_3d} steps)')
            m.make3D('mmff94', steps_3d)
            if steps_3d_localopt > 0:
                print(f'\tPerforming local optimization ({steps_3d_localopt} steps)')
                m.localopt('mmff94', steps_3d_localopt)

    # Prep for docking and write to file
    for idx, m in enumerate(mols):
        preparator = MoleculePreparation(hydrate=hydrate,
                                         keep_nonpolar_hydrogens=keep_nonpolar_hydrogens,
                                         pH_value=pH_value)
        preparator.prepare(m.OBMol)

        # preparator.write_pdbqt_file(f"{output_dir}/mol_{idx}.pdbqt")
        preparator.write_pdbqt_file(os.path.join(output_dir, f"mol_{idx}.pdbqt"))



# mol2
# sdf

# keep_nonpolar_hydrogens = False
# hydrate = False
# pH_value = None

# print(m)
#
# f'{m}'.split('\t')[0]

# Supported forcefields ['gaff', 'ghemical', 'mmff94', 'mmff94s', 'uff']

# obconversion = openbabel.OBConversion()
# obconversion.SetInFormat("sdf")
# obmol = openbabel.OBMol()
#
# notatend = obconversion.ReadFile(obmol, "example/RORytligand_testDerek.sdf")
# mol_nr = 1
# while notatend:
#     print(obmol.GetMolWt())
#     obmol = openbabel.OBMol()
#     notatend = obconversion.Read(obmol)
#
#     preparator = MoleculePreparation(hydrate=False)
#     preparator.prepare(obmol)
#     preparator.write_pdbqt_file(f'test_{mol_nr}.pdbqt')
#     mol_nr += 1
#
#
#


#     mols = [m for m in pybel.readfile(filename=infile, format=infile.split('.')[-1])]
#
#     for m in mols:
#         m.addh()
#
#     out = pybel.Outputfile(filename=outfile, format=outfile.split('.')[-1], overwrite=overwrite)
#     for mol in mols:
#         out.write(mol)
#     out.close()


# sanitize_mol2('1AZ8_lig.mol2', '1AZ8_lig_H.mol2')
# prepare_ligand('1AZ8_lig_H.mol2', '1AZ8_lig_H.pdbqt')

# 'example/RORytligand_testDerek.sdf'
#
# infile = '1AZ8_lig.mol2'
# mols = [m for m in pybel.readfile(filename='example/RORytligand_testDerek.sdf', format='sdf')]
#
# preparator = MoleculePreparation(hydrate=False)
# preparator.prepare(mols[0])
#
# preparator.write_pdbqt_file(outfile)
#
#
# sanitize_mol2('example/RORytligand_testDerek.sdf', 'example/RORytligand_testDerek.sdf', overwrite=True)
#
#
# mol = obutils.load_molecule_from_file('example/RORytligand_testDerek.sdf')
#
# preparator = MoleculePreparation(hydrate=False)
# preparator.prepare(mol)
#
# preparator.write_pdbqt_file(outfile)




