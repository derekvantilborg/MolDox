"""
Author: Derek van Tilborg
Data: dec 21 2021

Collection of functions to prepare a receptors pdb file for docking with Autodock Vina

"""


from openbabel import pybel, openbabel

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import MDAnalysis as MDA
from MDAnalysis.coordinates import PDB


def fix_protein(filename='', addHs_pH=7.4, output=None, renumber=False, keep_water=True, verbose=True):
    """ Function to clean/fix a pdb file and add hydrogens at a specific pH.
        Adapted from Angel J. Ruiz Moreno: https://github.com/AngelRuizMoreno/Jupyter_Dock under the MIT licence

    :param filename: (str) input .pdb file
    :param addHs_pH: (float) at which pH should we add hydrogens
    :param output: (str) name of output file
    :param renumber: (bool) if True, tries to renumber residues from 1 to n
    :param verbose: (bool) print out stuff if True
    :param keep_water: (bool) Keep waters if True

    """
    if output is None:
        output = f'{filename.split(".pdb")[0]}_H.pdb'
    fix = PDBFixer(filename=filename)

    # Cleanup pdb file
    if verbose:
        print(f"Cleanup pdb (replacing noncanonical residues and {'keeping' if keep_water else 'discarding'} waters")
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.removeHeterogens(keep_water)

    # Add missing atoms
    if verbose:
        print(f"Adding missing atoms if there are any")
    fix.findMissingAtoms()
    fix.addMissingAtoms()

    # Add missing hydrogens
    if verbose:
        print(f"Adding missing hydrogens (pH={addHs_pH})")
    fix.addMissingHydrogens(addHs_pH)

    # Write PDB file
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if renumber:
        try:
            if verbose:
                print('Renumbering residues')
            original = MDA.Universe(filename)
            from_fix = MDA.Universe(output)

            resNum = [res.resid for res in original.residues]
            for idx, res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save = PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Failed renumbering PDB residues')


def partial_charge(mol):
    """ Computate partial charges for a molecule. Uses GetPartialCharge to force computation of the Gasteiger
    partial charges.
    Adapted from Bharath Ramsundar: https://github.com/deepchem/deepchem under the MIT licence

    :param mol: (OBMol) Molecule which partial charges are computed

    """
    for atom in openbabel.OBMolAtomIter(mol):
        atom.GetPartialCharge()


def compute_partial_charges(input_pdb, output_pdbqt=None, verbose=False):
    """ Create a pdbqt file with partial charges
        Adapted from Bharath Ramsundar: https://github.com/deepchem/deepchem under the MIT licence

    :param input_pdb: (str) path to input pdb (hydrogenated)
    :param output_pdbqt: (str) name of output file (pdbqt file that has computed Gasteiger partial charges)
    :param verbose: (bool) prints out stuuf if True

    """
    # initiate OBConversion object
    charge_conversion = openbabel.OBConversion()

    if verbose:
        print("Convert input pdb to pdbqt.")
    charge_conversion.SetInAndOutFormats(str("pdb"), str("pdbqt"))

    if verbose:
        print("Make protein rigid.")
    charge_conversion.AddOption(str("c"), charge_conversion.OUTOPTIONS)
    charge_conversion.AddOption(str("r"), charge_conversion.OUTOPTIONS)

    if verbose:
        print("Preserve hydrogens")
    charge_conversion.AddOption(str("h"), charge_conversion.OUTOPTIONS)

    if verbose:
        print("Preserve atom indices")
    charge_conversion.AddOption(str("p"), charge_conversion.OUTOPTIONS)
    charge_conversion.AddOption(str("n"), charge_conversion.OUTOPTIONS)

    if verbose:
        print("Computing partial charges")
    mol = openbabel.OBMol()
    # read input file
    charge_conversion.ReadFile(mol, str(input_pdb))
    # compute partial charge
    partial_charge(mol)
    # write output file
    if output_pdbqt is None:
        output_pdbqt = f"{input_pdb}qt"
    charge_conversion.WriteFile(mol, str(output_pdbqt))


def prep_receptor(input_pdb, output_pdb=None, output_pdbqt=None, ph=7.4, keep_water=True, renumber=True, verbose=True):
    """ Function that takes care of fixing, hydrogenating, computing partial charge and conversion to pdbqt of a pdb"""

    fix_protein(filename=input_pdb,
                addHs_pH=ph,
                renumber=renumber,
                keep_water=keep_water,
                output=output_pdb,
                verbose=verbose)

    if output_pdb is None:
        output_pdb = f'{input_pdb.split(".pdb")[0]}_H.pdb'

    compute_partial_charges(input_pdb=output_pdb,
                            output_pdbqt=output_pdbqt,
                            verbose=verbose)
