
from docking.prepare import fix_protein
import prolif as plf
from prolif.plotting.network import LigNetwork
import MDAnalysis as mda


class InteractionMap:
    def __init__(self):
        pass

    def prep_protein(self):
        pass

    def mda(self):
        pass

    def load_ligands(self):
        pass

    def get_interaction_table(self):
        pass

    def get_interaction_network(self):
        pass

    def show(self):
        pass




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

