from meeko import MoleculePreparation
from meeko import obutils
from openbabel import pybel, openbabel
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


def sanitize_mol2(infile, outfile, overwrite=True):
    """

    :param infile: (str) pathname to input mol2 file
    :param outfile: (str) pathname to output mol2 file
    :param overwrite: (bool) if True, overwrite existing outfile
    """

    mols = [m for m in pybel.readfile(filename=infile, format=infile.split('.')[-1])]

    for m in mols:
        m.addh()

    out = pybel.Outputfile(filename=outfile, format=outfile.split('.')[-1], overwrite=overwrite)
    for mol in mols:
        out.write(mol)
    out.close()


def sanitize_protein():
    pass


def prepare_ligand(infile, outfile='1AZ8_lig_H.pdbqt'):
    mol = obutils.load_molecule_from_file(infile)

    preparator = MoleculePreparation(hydrate=False)
    preparator.prepare(mol)

    preparator.write_pdbqt_file(outfile)

def get_box():
    pass


def getbox(selection='sele', extending=6.0):

    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(selection)

    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)

    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX = (maxX + minX) / 2
    CenterY = (maxY + minY) / 2
    CenterZ = (maxZ + minZ) / 2

    cmd.delete('all')

    return {'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ}, \
           {'size_x': SizeX, 'size_y': SizeY, 'size_z': SizeZ}



def fix_protein(filename='', addHs_pH=7.4, output='', try_renumberResidues=False):
    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.removeHeterogens(True)
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if try_renumberResidues == True:
        try:
            original = mda.Universe(filename)
            from_fix = mda.Universe(output)

            resNum = [res.resid for res in original.residues]
            for idx, res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save = PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check excepton for extra details')


def force_partial_charge_computation(mol):
    """Force computation of partial charges for molecule.

    This function uses GetPartialCharge to force computation of the Gasteiger
    partial charges. This is an unfortunate hack, since it looks like the
    python openbabel API doesn't expose the OBGastChrg object which actually
    computes partial charges.

    Parameters
    ----------
    mol: OBMol
      Molecule on which we compute partial charges.
    """
    for obatom in openbabel.OBMolAtomIter(mol):
        obatom.GetPartialCharge()



def hydrogenate_and_compute_partial_charges(input_file, input_format,
                                            hyd_output=None,
                                            pdbqt_output=None,
                                            protein=True,
                                            verbose=False):
    """Outputs a hydrogenated pdb and a pdbqt with partial charges.

    Takes an input file in specified format. Generates two outputs:

    -) A pdb file that contains a hydrogenated (at pH 7.4) version of
       original compound.
    -) A pdbqt file that has computed Gasteiger partial charges. This pdbqt
       file is build from the hydrogenated pdb.

    Parameters
    ----------
    input_file: String
      Path to input file.
    input_format: String
      Name of input format.
    """

    # Since this function passes data to C++ obabel classes, we need to
    # constantly cast to str to convert unicode to char*
    if verbose:
        print("Create pdb with hydrogens added")
    hyd_conversion = openbabel.OBConversion()
    hyd_conversion.SetInAndOutFormats(str(input_format), str("pdb"))
    mol = openbabel.OBMol()
    hyd_conversion.ReadFile(mol, str(input_file))
    # AddHydrogens(not-polaronly, correctForPH, pH)
    mol.AddHydrogens(False, True, 7.4)
    hyd_conversion.WriteFile(mol, str(hyd_output))

    if verbose:
        print("Create a pdbqt file from the hydrogenated pdb above.")
    charge_conversion = openbabel.OBConversion()
    charge_conversion.SetInAndOutFormats(str("pdb"), str("pdbqt"))

    if protein and verbose:
        print("Make protein rigid.")
    if protein:
        charge_conversion.AddOption(str("c"), charge_conversion.OUTOPTIONS)
        charge_conversion.AddOption(str("r"), charge_conversion.OUTOPTIONS)
    if verbose:
        print("Preserve hydrogens")
    charge_conversion.AddOption(str("h"), charge_conversion.OUTOPTIONS)
    if verbose:
        print("Preserve atom indices")
    charge_conversion.AddOption(str("p"), charge_conversion.OUTOPTIONS)
    if verbose:
        print("preserve atom indices.")
    charge_conversion.AddOption(str("n"), charge_conversion.OUTOPTIONS)

    if verbose:
        print("About to run obabel conversion.")
    mol = openbabel.OBMol()
    charge_conversion.ReadFile(mol, str(hyd_output))
    force_partial_charge_computation(mol)
    charge_conversion.WriteFile(mol, str(pdbqt_output))


def pdbqt_to_sdf(pdbqt_file=None, output=None):

    results = [m for m in pybel.readfile(filename=pdbqt_file, format='pdbqt')]
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)
    for pose in results:

        pose.data.update({'Pose': pose.data['MODEL']})
        pose.data.update({'Score': pose.data['REMARK'].split()[2]})
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']

        out.write(pose)
    out.close()
