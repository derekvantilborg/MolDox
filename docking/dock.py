from rdkit import Chem
from pymol import cmd
from vina import Vina
from docking.prepare_receptor import prep_receptor
from docking.prepare_ligands import hydrogenate_mol_from_file, prep_mol_from_file, prep_ligands
from docking.utils import find_box, pdbqt_to_sdf
from docking.viewer import Mol3D, InteractionMap
import os
from os.path import basename


class MolDox:
    def __init__(self, output_dir, receptor=None, ref_ligand=None):

        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        self.receptor = receptor
        self.receptor_pdbqt = None
        self.receptor_hydrogenated = None

        self.ref_ligand = ref_ligand
        self.ref_ligand_hydrogenated = None
        self.ref_ligand_pdbqt = None

        self.ligands = None
        self.ligands_pdbqt = {}

        self.box_center = None
        self.box_size = None

        self.docking_results = {}

        self.view = None

    def prep_reference_ligand(self, ref_ligand=None, ref_ligand_pdbqt=None, ref_ligand_hydrogenated=None):

        # Set and remember the pathnames of input/intermediate/output files
        if self.ref_ligand is None:
            self.ref_ligand = ref_ligand

        if ref_ligand_hydrogenated is None:
            # self.ref_ligand_hydrogenated = f'{self.output_dir}/{basename(self.ref_ligand).split(".")[0]}_H.mol2'
            self.ref_ligand_hydrogenated = os.path.join(self.output_dir, f'{basename(self.ref_ligand).split(".")[0]}_H.mol2')
        else:
            # self.ref_ligand_hydrogenated = f'{self.output_dir}/{ref_ligand_hydrogenated}'
            self.ref_ligand_hydrogenated = os.path.join(self.output_dir, ref_ligand_hydrogenated)

        if ref_ligand_pdbqt is None:
            self.ref_ligand_pdbqt = f"{self.ref_ligand_hydrogenated}".replace('mol2', 'pdbqt')
        else:
            # self.ref_ligand_pdbqt = f'{self.output_dir}/{ref_ligand_pdbqt}'
            self.ref_ligand_pdbqt = os.path.join(self.output_dir, ref_ligand_pdbqt)

        # test_output/1AZ8_lig.mol2_H.pdb

        hydrogenate_mol_from_file(self.ref_ligand, self.ref_ligand_hydrogenated)
        prep_mol_from_file(self.ref_ligand_hydrogenated, self.ref_ligand_pdbqt)

    def prep_ligands(self, infile, format='sdf', hydrogenate=True, hydrate=False, keep_nonpolar_hydrogens=False,
                     pH_value=None, make3d=True, remake3d=False, forcefield='mmff94', steps_3d=50,
                     steps_3d_localopt=500):

        self.ligands = infile

        # if not os.path.exists(f"{self.output_dir}/ligands"):
        #     os.mkdir(f"{self.output_dir}/ligands")
        if not os.path.exists(os.path.join(self.output_dir, 'ligands')):
            os.mkdir(os.path.join(self.output_dir, 'ligands'))



        prep_ligands(infile,
                     format=format,
                     # output_dir=f"{self.output_dir}/ligands",
                     output_dir=os.path.join(self.output_dir, 'ligands'),
                     hydrogenate=hydrogenate,
                     hydrate=hydrate,
                     keep_nonpolar_hydrogens=keep_nonpolar_hydrogens,
                     pH_value=pH_value,
                     make3d=make3d,
                     remake3d=remake3d,
                     forcefield=forcefield,
                     steps_3d=steps_3d,
                     steps_3d_localopt=steps_3d_localopt)

        # Add all pdbqt files of the prepped ligands to a dict
        # for pdbqt in os.listdir(f"{self.output_dir}/ligands"):
        #     self.ligands_pdbqt[pdbqt.split('.')[0]] = f"{self.output_dir}/ligands/{pdbqt}"
        for pdbqt in os.listdir(os.path.join(self.output_dir, 'ligands')):
            self.ligands_pdbqt[pdbqt.split('.')[0]] = os.path.join(self.output_dir, 'ligands', pdbqt)

    def prep_receptor(self, receptor=None, ligand=None, receptor_hydrogenated=None, receptor_pdbqt=None,
                      ph=7.4, keep_water=True, renumber=True, verbose=True):

        # Set and remember the pathnames of input/intermediate/output files
        if self.receptor is None:
            self.receptor = receptor

        if receptor_hydrogenated is None:
            # self.receptor_hydrogenated = f'{self.output_dir}/{basename(self.receptor).split(".pdb")[0]}_H.pdb'
            self.receptor_hydrogenated = os.path.join(self.output_dir, f'{basename(self.receptor).split(".pdb")[0]}_H.pdb')
        else:
            # self.receptor_hydrogenated = f'{self.output_dir}/{receptor_hydrogenated}'
            self.receptor_hydrogenated = os.path.join(self.output_dir, receptor_hydrogenated)

        if receptor_pdbqt is None:
            # self.receptor_pdbqt = f"{self.output_dir}/{basename(self.receptor_hydrogenated)}qt"
            self.receptor_pdbqt = os.path.join(self.output_dir, f"{basename(self.receptor_hydrogenated)}qt")
        else:
            # self.receptor_pdbqt = f'{self.output_dir}/{receptor_pdbqt}'
            self.receptor_pdbqt = os.path.join(self.output_dir, receptor_pdbqt)

        # Actually clean the receptor pdb
        prep_receptor(input_pdb=self.receptor,
                      output_pdb=self.receptor_hydrogenated,
                      output_pdbqt=self.receptor_pdbqt,
                      ph=ph,
                      keep_water=keep_water,
                      renumber=renumber,
                      verbose=verbose)

        if ligand is not None:
            self.prep_reference_ligand(ligand)

    def find_box(self, receptor=None, ligand=None, receptor_format='pdb', ligand_format='mol2', box_extension=5):

        if receptor is None:
            receptor = self.receptor_hydrogenated
        if ligand is None:
            ligand = self.ref_ligand_hydrogenated


        self.box_center, self.box_size = find_box(receptor=receptor,
                                                  ligand=ligand,
                                                  receptor_format=receptor_format,
                                                  ligand_format=ligand_format,
                                                  box_extension=box_extension)

    def autodock_all(self, receptor_pdbqt=None, box_center=None, box_size=None, box_extension=5,
                     exhaustiveness=10, n_poses=10, scoring_function='vina', cpu_cores=0, seed=42):

        if box_center is None or box_size is None:
            if self.box_center is None or self.box_size is None:
                self.find_box(box_extension=box_extension)
            box_center = [self.box_center['center_x'], self.box_center['center_y'], self.box_center['center_z']]
            box_size = [self.box_size['size_x'], self.box_size['size_y'], self.box_size['size_z']]

        if receptor_pdbqt is None:
            receptor_pdbqt = self.receptor_pdbqt

        # if not os.path.exists(f"{self.output_dir}/docking_results"):
        #     os.mkdir(f"{self.output_dir}/docking_results")

        if not os.path.exists(os.path.join(self.output_dir, 'docking_results')):
            os.mkdir(os.path.join(self.output_dir, 'docking_results'))



        # Perform the docking for every ligand
        print(f"Docking {len(self.ligands_pdbqt)} ligands")
        for lig, lig_pdbqt in self.ligands_pdbqt.items():
            autodock_vina(receptor_pdbqt=receptor_pdbqt,
                          ligand_pdbqt=lig_pdbqt,
                          # output_file=f"{self.output_dir}/docking_results/{lig}.sdf",
                          output_file=os.path.join(self.output_dir, 'docking_results', f"{lig}.sdf"),
                          box_center=box_center,
                          box_size=box_size,
                          exhaustiveness=exhaustiveness,
                          n_poses=n_poses,
                          scoring_function=scoring_function,
                          cpu_cores=cpu_cores,
                          seed=seed)

            # self.docking_results[lig] = f"{self.output_dir}/docking_results/{lig}.sdf"
            self.docking_results[lig] = os.path.join(self.output_dir, 'docking_results', f"{lig}.sdf")


    def dock(self, ligand_pdbqt=None, receptor_pdbqt=None, output_file=None, box_center=None, box_size=None,
             exhaustiveness=10, n_poses=10, scoring_function='vina', cpu_cores=0, seed=42):

        # If no ligand is provided, redock the reference ligand
        if ligand_pdbqt is None:
            ligand_pdbqt = self.ref_ligand_pdbqt

        if receptor_pdbqt is None:
            receptor_pdbqt = self.receptor_pdbqt

        if box_center is None:
            box_center = [self.box_center['center_x'], self.box_center['center_y'], self.box_center['center_z']]
        if box_size is None:
            box_size = [self.box_size['size_x'], self.box_size['size_y'], self.box_size['size_z']]

        # If no output name is given, use a generic name
        if output_file is None:
            # output_file = f"{self.output_dir}/docking_results.sdf"
            output_file = os.path.join(self.output_dir, "docking_results.sdf")

        # Perform the docking
        autodock_vina(receptor_pdbqt, ligand_pdbqt, output_file, box_center, box_size,
                      exhaustiveness=exhaustiveness, n_poses=n_poses,
                      scoring_function=scoring_function, cpu_cores=cpu_cores, seed=seed)

    def interactions(self, docking_results_sdf):

        if not docking_results_sdf.endswith('.sdf'):
            docking_results_sdf = self.docking_results[docking_results_sdf]

        self.interaction_map = InteractionMap(self.receptor_hydrogenated, docking_results_sdf)
        return self.interaction_map.show()

    def viewer(self, receptor=None, ligand=None, receptor_color='white', surface_opacity=0.75, ligand_color='cyan',
               ligand_colorscheme=None, ligand_opacity=1):

        if receptor is None:
            receptor = self.receptor_hydrogenated
        if ligand is None:
            ligand = self.ref_ligand_hydrogenated

        self.view = Mol3D()
        self.view.add_layer(receptor)
        self.view.add_surface(color=receptor_color, opacity=surface_opacity)
        self.view.cartoon()
        self.view.add_layer(ligand)
        self.view.sticks(color=ligand_color, colorscheme=ligand_colorscheme, opacity=ligand_opacity)
        self.view.add_outline()

    def viewer_add_ligand(self, ligand, pose=0, color='orange', colorscheme=None, opacity=1):
        if self.view is None:
            self.viewer()

        if not ligand.endswith('.sdf'):
            ligand = self.docking_results[ligand]

        results=Chem.SDMolSupplier(ligand)
        p = Chem.MolToMolBlock(results[pose], False)

        self.view.add_layer_from_object(p)
        self.view.sticks(color=color, colorscheme=colorscheme, opacity=opacity)

    def show(self):
        if self.view is None:
            self.viewer()
        self.view.show()


def autodock_vina(receptor_pdbqt, ligand_pdbqt, output_file, box_center, box_size, exhaustiveness=10, n_poses=10,
                  scoring_function='vina', cpu_cores=0, seed=42):
    """

    The vina scoring function (lower is better) is a hybrid function described in doi: 10.1002/jcc.21334, PMID: 19499576


    :param receptor_pdbqt: (str) path to receptor pdbqt file
    :param ligand_pdbqt: (str) path to ligand pdbqt file
    :param output_file: (str) path of output file (func will create a .pdbqt and .sdf file)
    :param box_center: (list) [center_x, center_y, center_z]
    :param box_size: (list) [size_x, size_y, size_z]
    :param exhaustiveness: (int) Number of MC run (default: 10)
    :param n_poses: (int) number of pose to generate (default: 20)
    :param scoring_function: (str) Scoring method to score docked ligands: vina, ad4 (default: vina)
    :param cpu_cores: (int) number of cpu cores to use (default=0, which means all cores)
    :param seed: (int) random seed

    """

    v = Vina(sf_name=scoring_function, cpu=cpu_cores, seed=seed)

    # Set the receptor and ligand
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    # Define the docking maps using the box center and size
    v.compute_vina_maps(center=box_center, box_size=box_size)

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
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(f"{output_file.split('.')[0]}.pdbqt", n_poses=n_poses, overwrite=True)

    # Docking results to SDF #
    pdbqt_to_sdf(pdbqt_file=f"{output_file.split('.')[0]}.pdbqt", output=f"{output_file.split('.')[0]}.sdf")


#
# v.compute_vina_maps(center=[self.box_center['center_x'],
#                             self.box_center['center_y'],
#                             self.box_center['center_z']],
#                     box_size=[self.box_size['size_x'],
#                               self.box_size['size_y'],
#                               self.box_size['size_z']])