def move_some_atoms(dfc, atoms_idcs_list, dx=0, dy=0, dz=0):
    for atom in dfc.atoms:
        if atom['atom_id'] in atoms_idcs_list:
            atom['x'] += dx
            atom['y'] += dy
            atom['z'] += dz
