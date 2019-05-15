from datafile_content import DatafileContent
from datafile_content_modifiers import center_in_box, get_bbox


def datafile_content_extract(dfc, atom_ids_list):
    '''
    Removes everything from the dfc (DatafileContent class) except
    atoms having ids that are present in atom_ids_list (list type).
    '''
    atom_idcs_to_remove = []
    for atom_idx, atom in enumerate(dfc.atoms):
        if atom['atom_id'] not in atom_ids_list:
            atom_idcs_to_remove.append(atom_idx)
    for atom_idx in sorted(atom_idcs_to_remove)[::-1]:
        dfc.atoms.pop(atom_idx)
    print('removed {0} atoms, left {1}'.format(
        len(atom_idcs_to_remove), len(dfc.atoms)))

    bond_idcs_to_remove = []
    for bond_idx, bond in enumerate(dfc.bonds):
        if any([bond['atom_one_id'] not in atom_ids_list,
                bond['atom_two_id'] not in atom_ids_list]):
            bond_idcs_to_remove.append(bond_idx)
    for bond_idx in sorted(bond_idcs_to_remove)[::-1]:
        dfc.bonds.pop(bond_idx)
    print('removed {0} bonds, left {1}'.format(
        len(bond_idcs_to_remove), len(dfc.bonds)))

    angle_idcs_to_remove = []
    for angle_idx, angle in enumerate(dfc.angles):
        if any([angle['atom_one_id'] not in atom_ids_list,
                angle['atom_two_id'] not in atom_ids_list,
                angle['atom_three_id'] not in atom_ids_list]):
            angle_idcs_to_remove.append(angle_idx)
    for angle_idx in sorted(angle_idcs_to_remove)[::-1]:
        dfc.angles.pop(angle_idx)
    print('removed {0} angles, left {1}'.format(
        len(angle_idcs_to_remove), len(dfc.angles)))

    dihedral_idcs_to_remove = []
    for dihedral_idx, dihedral in enumerate(dfc.dihedrals):
        if any([dihedral['atom_one_id'] not in atom_ids_list,
                dihedral['atom_two_id'] not in atom_ids_list,
                dihedral['atom_three_id'] not in atom_ids_list,
                dihedral['atom_four_id'] not in atom_ids_list]):
            dihedral_idcs_to_remove.append(dihedral_idx)
    for dihedral_idx in sorted(dihedral_idcs_to_remove)[::-1]:
        dfc.dihedrals.pop(dihedral_idx)
    print('removed {0} dihedrals, left {1}'.format(
        len(dihedral_idcs_to_remove), len(dfc.dihedrals)))

    improper_idcs_to_remove = []
    for improper_idx, improper in enumerate(dfc.impropers):
        if any([improper['atom_one_id'] not in atom_ids_list,
                improper['atom_two_id'] not in atom_ids_list,
                improper['atom_three_id'] not in atom_ids_list,
                improper['atom_four_id'] not in atom_ids_list]):
            improper_idcs_to_remove.append(improper_idx)
    for improper_idx in sorted(improper_idcs_to_remove)[::-1]:
        dfc.impropers.pop(improper_idx)
    print('removed {0} impropers, left {1}'.format(
        len(improper_idcs_to_remove), len(dfc.impropers)))
    dfc.atoms_count = len(dfc.atoms)
    dfc.bonds_count = len(dfc.bonds)
    dfc.angles_count = len(dfc.angles)
    dfc.dihedrals_count = len(dfc.dihedrals)
    dfc.impropers_count = len(dfc.impropers)
    if dfc.atoms_count == 0:
        dfc.atoms = None
        dfc.atom_types = None
        dfc.bonds = None
        dfc.bond_types = None
        dfc.angles = None
        dfc.angle_types = None
        dfc.dihedrals = None
        dfc.dihedral_types = None
        dfc.impropers = None
        dfc.improper_types = None
        print('Warning! no atoms left after datafile_content_extract call')
    if dfc.bonds_count == 0:
        dfc.bonds = None
        dfc.bond_types = None
    if dfc.angles_count == 0:
        dfc.angles = None
        dfc.angle_types = None
    if dfc.dihedrals_count == 0:
        dfc.dihedrals = None
        dfc.dihedral_types = None
    if dfc.impropers_count == 0:
        dfc.impropers = None
        dfc.improper_types = None


if __name__ == '__main__':
    fname = 'mix.data'
    dfc = DatafileContent(fname)

    #atom_ids_list = range(1561, 1561+192)  # leave single modifier molecule
    atom_ids_list = range(1, 721)  # leave only mmt

    datafile_content_extract(dfc, atom_ids_list)

    dfc.reassign_atom_ids()
    dfc.reassign_bond_ids()
    dfc.reassign_angle_ids()
    dfc.reassign_dihedral_ids()
    dfc.reassign_improper_ids()

    dfc.reassign_atom_ids()
    dfc.reassign_atom_types()
    dfc.reassign_bond_types()
    dfc.reassign_angle_types()
    dfc.reassign_dihedral_types()

    center_in_box(dfc, border=1)
    dfc.zlo -= 3

    dfc.write('mmt.data')
