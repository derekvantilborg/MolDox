from rdkit import Chem
from pymol import cmd
from vina import Vina
from docking.prepare_receptor import prep_receptor
from docking.prepare_ligands import prep_single_ligand
from docking.utils import find_box, pdbqt_to_sdf
from viewer.mol3d import Viewer
import os


class MolDox:
    def __init__(self, output_dir, receptor=None, ligand=None):

        self.output_dir = output_dir

        self.receptor = receptor
        self.receptor_pdbqt = None
        self.receptor_hydrogenated = None

        self.ligand = ligand
        self.ligands_hydrogenated = {}
        self.ligands_pdbqt = {}

        self.docking_results = {}





    def prep_ligands(self):
        pass

    def prep_receptor(self, receptor=None, output_pdb=None, output_pdbqt=None, ph=7.4, keep_water=True, renumber=True,
                      verbose=True):

        # Set and remember the pathnames of input/intermediate/output files
        if self.receptor is None:
            self.receptor = receptor

        if output_pdb is not None:
            self.receptor_hydrogenated = f'{self.receptor.split(".pdb")[0]}_H.pdb'
        else:
            self.receptor_hydrogenated = output_pdb

        if output_pdbqt is not None:
            self.receptor_pdbqt = f"{self.receptor_hydrogenated}qt"
        else:
            self.receptor_pdbqt = output_pdb

        # Actually clean the receptor pdb
        prep_receptor(input_pdb=self.receptor,
                      output_pdb=self.receptor_hydrogenated,
                      output_pdbqt=self.receptor_pdbqt,
                      ph=ph,
                      keep_water=keep_water,
                      renumber=renumber,
                      verbose=verbose)

    def find_box(self):
        pass
        # self.box_center, self.box_size = find_box(receptor=self.cleaned_receptor,
        #                                           ligand=f'{output_dir}/1AZ8_lig_H.mol2',
        #                                           receptor_format=self.cleaned_receptor.split('.')[-1],
        #                                           ligand_format='mol2', box_extension=5)

    def dock(self):
        pass

    def view(self, ligand=None):
        pass
