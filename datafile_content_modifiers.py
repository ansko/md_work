from datafile_content import DatafileContent


def translate_in_box(dfc, dx=None, dy=None, dz=None):
    '''
    Translate all atoms without moving the box 
    '''
    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    if len(set([dx, dy, dz, None])) == 1:
        return
    for atom in dfc.atoms:
        if dx is not None:
            atom['x'] += dx
        if dy is not None:
            atom['y'] += dy
        if dz is not None:
            atom['z'] += dz
    print('translating in box...')

def get_bbox(dfc):
    '''
    Return bounding box as a dict
    '''
    bbox_xlo = bbox_xhi = dfc.atoms[0]['x']
    bbox_ylo = bbox_yhi = dfc.atoms[0]['y']
    bbox_zlo = bbox_zhi = dfc.atoms[0]['z']
    for atom in dfc.atoms[1:]:
        bbox_xlo = min(bbox_xlo, atom['x'])
        bbox_xhi = max(bbox_xhi, atom['x'])
        bbox_ylo = min(bbox_ylo, atom['y'])
        bbox_yhi = max(bbox_yhi, atom['y'])
        bbox_zlo = min(bbox_zlo, atom['z'])
        bbox_zhi = max(bbox_zhi, atom['z'])
    return {
        'bbox_xlo': bbox_xlo, 'bbox_xhi': bbox_xhi,
        'bbox_ylo': bbox_ylo, 'bbox_yhi': bbox_yhi,
        'bbox_zlo': bbox_zlo, 'bbox_zhi': bbox_zhi}

def center_in_box(dfc, border=0):
    '''
    Make all atom flags (nx, ny, nz) equal to zero.
    border is additional emptiness around atom (to avoid close
    contacts when applying periodic boundary conditions).
    '''
    lx = dfc.xhi - dfc.xlo
    ly = dfc.yhi - dfc.ylo
    lz = dfc.zhi - dfc.zlo
    for atom in dfc.atoms:
        atom['x'] += atom['nx'] * lx
        atom['y'] += atom['ny'] * ly
        atom['z'] += atom['nz'] * lz
        atom['nx'] = atom['ny'] = atom['nz'] = 0
    bbox = get_bbox(dfc)
    dfc.xlo = bbox['bbox_xlo'] - border; dfc.xhi = bbox['bbox_xhi'] + border
    dfc.ylo = bbox['bbox_ylo'] - border; dfc.yhi = bbox['bbox_yhi'] + border
    dfc.zlo = bbox['bbox_zlo'] - border; dfc.zhi = bbox['bbox_zhi'] + border 


if __name__ == '__main__':
    fname = 'data_structures/modifier.data'
    dfc = DatafileContent(fname)

    dfc.reassign_atom_ids()
    dfc.reassign_atom_types()
    dfc.reassign_bond_types()
    dfc.reassign_angle_types()
    dfc.reassign_dihedral_types()

    #translate_in_box(dfc, dy=10)
    center_in_box(dfc)

    bbox = get_bbox(dfc)
    dfc.xlo = bbox['bbox_xlo'] - 1; dfc.xhi = bbox['bbox_xhi'] + 1
    dfc.ylo = bbox['bbox_ylo'] - 1; dfc.yhi = bbox['bbox_yhi'] + 1
    dfc.zlo = bbox['bbox_zlo'] - 1; dfc.zhi = bbox['bbox_zhi'] + 1

    dfc.write('new_mod.data')
