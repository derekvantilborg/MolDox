
from meeko import MoleculePreparation
from meeko import obutils
from openbabel import pybel, openbabel
import os


def sanitize_mol(infile, outfile, overwrite=True):
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


def prepare_ligand(infile, output_pdbqt):
    mol = obutils.load_molecule_from_file(infile)

    preparator = MoleculePreparation(hydrate=False)
    preparator.prepare(mol)

    preparator.write_pdbqt_file(output_pdbqt)


def prep_single_ligand(infile, output_pdbqt):
    """ Sanitize and hydrate a single molecule from a file. Create a pdbqt file suitable for Autodock Vina"""

    sanitize_mol(infile, f"{infile.split('.')[-2]}_H.mol2")
    prepare_ligand(f"{infile.split('.')[-2]}_H.mol2", output_pdbqt)


def prep_ligands(infile, output_dir='.'):

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    mols = [m for m in pybel.readfile(filename=infile, format=infile.split('.')[-1])]

    # Add hydrogens
    for m in mols:
        m.addh()

    for idx, mol in enumerate(mols):

        # Write to output file
        out = pybel.Outputfile(filename=f"{output_dir}/{infile.split('.')[-2]}_temp_H.mol2",
                               format='mol2', overwrite=True)
        out.write(mol)
        out.close()

        prepare_ligand(f"{output_dir}/{infile.split('.')[-2]}_temp_H.mol2",
                       f"{output_dir}/{infile.split('.')[-2]}_H_{idx}.pdbqt")

        os.remove(f"{output_dir}/{infile.split('.')[-2]}_temp_H.mol2")


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




