
from openbabel import pybel
from pymol import cmd


def find_box(receptor, ligand, receptor_format='pdb', ligand_format='mol2', box_extension=5):
    """ Get the box coordinates from the ligand and pdb

    {'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ}, \
           {'size_x': SizeX, 'size_y': SizeY, 'size_z': SizeZ}

    :param receptor: (str) filename of receptor
    :param ligand: (str) filename of ligand
    :param receptor_format: (str) Filetype of receptor (default=pdb)
    :param ligand_format: (str) Filetype of ligand (default=mol2)
    :param box_extension: (int/float) Extend box of ligand by n Angstrom (default=5)

    :return: {'center_x': float, 'center_y': float, 'center_z': float},
             {'size_x': float, 'size_y': float, 'size_z': float}
    """

    cmd.load(filename=receptor, format=receptor_format, object='prot')
    cmd.load(filename=ligand, format=ligand_format, object='lig')

    center, size = getbox(selection='lig', extending=box_extension)

    cmd.delete('all')

    return center, size


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


def pdbqt_to_sdf(pdbqt_file=None, output=None):

    results = [m for m in pybel.readfile(filename=pdbqt_file, format='pdbqt')]
    out = pybel.Outputfile(filename=output, format='sdf', overwrite=True)
    for pose in results:

        pose.data.update({'Pose': pose.data['MODEL']})
        pose.data.update({'Score': pose.data['REMARK'].split()[2]})
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']

        out.write(pose)
    out.close()
