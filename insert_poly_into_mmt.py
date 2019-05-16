import copy

from datafile_content import DatafileContent
from datafile_content_modifiers import get_bbox


mmt_thickness = 8


def extend_mmt_util(dfc_mmt, dz=0):
    '''
    Extend the cell (containing only mmt) in z direction
    '''
    assert(dz >= 0)  # may be changed later
    zlo = dfc_mmt.zlo
    zhi = dfc_mmt.zhi
    atoms_zlo = dfc_mmt.atoms[0]['z']
    atoms_zhi = dfc_mmt.atoms[0]['z']
    for atom in dfc_mmt.atoms[1:]:
        atoms_zlo = min(atoms_zlo, atom['z'])
        atoms_zhi = max(atoms_zhi, atom['z'])
    if atoms_zhi - atoms_zlo < mmt_thickness:  # mmt does not cross box along z
        dz_lo = abs(dfc_mmt.atoms[0]['z'] - dfc_mmt.zlo)
        dz_hi = abs(dfc_mmt.atoms[0]['z'] - dfc_mmt.zhi)
        if dz_lo < dz_hi:  # mmt is near low border -> extend high border
            dfc_mmt.zhi += dz
        else:
            dfc_mmt.zlo -= dz


def extend_mmt_ang_mod_util(dfc_mmt, dz=0):
    '''
    Extend the cell (containing bot mmt and modifier) in z direction
    
    '''


    # TODO : everything is not ready!
    assert(dz >= 0)  # may be changed later
    zlo = dfc_mmt.zlo
    zhi = dfc_mmt.zhi
    atoms_zlo = dfc_mmt.atoms[0]['z']
    atoms_zhi = dfc_mmt.atoms[0]['z']
    for atom in dfc_mmt.atoms[1:]:
        atoms_zlo = min(atoms_zlo, atom['z'])
        atoms_zhi = max(atoms_zhi, atom['z'])
    if atoms_zhi - atoms_zlo < mmt_thickness:  # mmt does not cross box along z
        dz_lo = abs(dfc_mmt.atoms[0]['z'] - dfc_mmt.zlo)
        dz_hi = abs(dfc_mmt.atoms[0]['z'] - dfc_mmt.zhi)
        if dz_lo < dz_hi:  # mmt is near low border -> extend high border
            dfc_mmt.zhi += dz
        else:
            dfc_mmt.zlo -= dz


def util_get_poly_ns(dfc_mmt, dfc_poly, border=0):
    bbox_mmt = get_bbox(dfc_mmt)
    bbox_poly = get_bbox(dfc_poly)
    mmt_lx = dfc_mmt.xhi - dfc_mmt.xlo
    mmt_ly = dfc_mmt.yhi - dfc_mmt.ylo
    mmt_lz = dfc_mmt.zhi - dfc_mmt.zlo
    poly_cell_lx = dfc_poly.xhi - dfc_poly.xlo
    poly_cell_ly = dfc_poly.yhi - dfc_poly.ylo
    poly_cell_lz = dfc_poly.zhi - dfc_poly.zlo
    poly_bbox_lx = bbox_poly['bbox_xhi'] - bbox_poly['bbox_xlo']
    poly_bbox_ly = bbox_poly['bbox_yhi'] - bbox_poly['bbox_ylo']
    poly_bbox_lz = bbox_poly['bbox_zhi'] - bbox_poly['bbox_zlo']
    dxlo = bbox_mmt['bbox_xlo'] - bbox_poly['bbox_xlo']
    dylo = bbox_mmt['bbox_ylo'] - bbox_poly['bbox_ylo']
    dzlo = bbox_mmt['bbox_zlo'] - bbox_poly['bbox_zlo']
    dxhi = bbox_mmt['bbox_xhi'] - bbox_poly['bbox_xhi']
    dyhi = bbox_mmt['bbox_yhi'] - bbox_poly['bbox_yhi']
    dzhi = bbox_mmt['bbox_zhi'] - bbox_poly['bbox_zhi']
    mmt_cell_lo_bbox_lo = dfc_mmt.zlo - bbox_mmt['bbox_zlo']
    mmt_cell_hi_bbox_hi = bbox_mmt['bbox_zhi'] - dfc_mmt.zhi
    bbox_mmt_empty = {
        'bbox_xlo': bbox_mmt['bbox_xlo'], 'bbox_xhi': bbox_mmt['bbox_xhi'],
        'bbox_ylo': bbox_mmt['bbox_ylo'], 'bbox_yhi': bbox_mmt['bbox_yhi'],
    }
    if mmt_cell_lo_bbox_lo > mmt_cell_hi_bbox_hi:
        bbox_mmt_empty['bbox_zlo'] = bbox_mmt['bbox_zhi']
        bbox_mmt_empty['bbox_zhi'] = dfc_mmt.zhi
    else:
        bbox_mmt_empty['bbox_zlo'] = dfc_mmt.zlo
        bbox_mmt_empty['bbox_zhi'] = bbox_mmt['bbox_zlo']
    if mmt_lx < poly_bbox_lx:
        print('error! small box in x direction!')
        print('mmt lx:', mmt_lx)
        print('poly bbox lx:', poly_bbox_lx)
        return False
    if mmt_ly < poly_bbox_ly:
        print('error! small box in y direction!')
        print('mmt ly', mmt_ly)
        print('poly bbox ly:', poly_bbox_ly)
        return False
    if mmt_lz < poly_bbox_lz:
        print('error! small box in z direction!')
        print('mmt lz:', mmt_lz)
        print('poly bbox lz:', poly_bbox_lz)
        return False
    empty_lz = bbox_mmt_empty['bbox_zhi'] - bbox_mmt_empty['bbox_zlo']
    # how much polys can be inserted along axes
    #print(empty_lz, poly_bbox_lz)
    n_poly_x = int(mmt_lx // (poly_bbox_lx + 2*border))
    n_poly_y = int(mmt_ly // (poly_bbox_ly + 2*border))
    n_poly_z = int(empty_lz // (poly_bbox_lz + 2*border))
    return n_poly_x, n_poly_y, n_poly_z


def insert_poly_into_mmt(dfc_mmt, dfc_poly, border=0):
    '''
    Insert (tryes to insert) polymer into empty mmt interlayer space
    return True if suceed otherwise False
    border is additional emptiness between poly's bbox-es
    '''
    def util_compare_atom_types(old, new, float_compare_accuracy=1e-5):
        '''
        return True if are the same (except comments maybe) otherwise False
        '''
        d_mass =  old['mass'] - new['mass']
        s_mass =  old['mass'] + new['mass']
        if abs(d_mass / s_mass) > float_compare_accuracy:
            return False
        d_eps = old['pair_coeffs']['coeff_eps'] - new['pair_coeffs']['coeff_eps']
        s_eps = old['pair_coeffs']['coeff_eps'] + new['pair_coeffs']['coeff_eps']
        if d_eps == 0 and s_eps == 0:
            ...
        elif abs(d_eps / s_eps) > float_compare_accuracy:
            return False
        d_sig = old['pair_coeffs']['coeff_sig'] - new['pair_coeffs']['coeff_sig']
        s_sig = old['pair_coeffs']['coeff_sig'] + new['pair_coeffs']['coeff_sig']
        if d_sig == 0 and s_sig == 0:
            ...
        elif abs(d_sig / s_sig) > float_compare_accuracy:
            return False
        old_label = old['pair_coeffs']['comment']
        new_label = new['pair_coeffs']['comment']
        if ((old_label and new_label) and (old_label != new_label)):
            return False
        return True

    def util_compare_bond_types(old, new, float_compare_accuracy=1e-5):
        '''
        return True if are the same (except comments maybe) otherwise False
        '''
        d_k = old['bond_coeffs']['coeff_k'] - new['bond_coeffs']['coeff_k']
        s_k = old['bond_coeffs']['coeff_k'] + new['bond_coeffs']['coeff_k']
        if d_k == 0 and s_k == 0:
            ...
        elif abs(d_k / s_k) > float_compare_accuracy:
            return False
        d_l = old['bond_coeffs']['coeff_l'] - new['bond_coeffs']['coeff_l']
        s_l = old['bond_coeffs']['coeff_l'] + new['bond_coeffs']['coeff_l']
        if d_l == 0 and s_l == 0:
            ...
        elif abs(d_l / s_l) > float_compare_accuracy:
            return False
        old_label = old['bond_coeffs']['comment']
        new_label = new['bond_coeffs']['comment']
        if ((old_label and new_label) and (old_label != new_label)):
            return False
        return True

    def util_compare_angle_types(old, new, float_compare_accuracy=1e-5):
        '''
        return True if are the same (except comments maybe) otherwise False
        '''
        d_k = old['angle_coeffs']['coeff_k'] - new['angle_coeffs']['coeff_k']
        s_k = old['angle_coeffs']['coeff_k'] + new['angle_coeffs']['coeff_k']
        if d_k == 0 and s_k == 0:
            ...
        elif abs(d_k / s_k) > float_compare_accuracy:
            return False
        d_theta = (old['angle_coeffs']['coeff_theta']
               - new['angle_coeffs']['coeff_theta'])
        s_theta = (old['angle_coeffs']['coeff_theta']
               + new['angle_coeffs']['coeff_theta'])
        if d_theta == 0 and s_theta == 0:
            ...
        elif abs(d_theta / s_theta) > float_compare_accuracy:
            return False
        old_label = old['angle_coeffs']['comment']
        new_label = new['angle_coeffs']['comment']
        if ((old_label and new_label) and (old_label != new_label)):
            return False
        return True

    def util_compare_dihedral_types(old, new, float_compare_accuracy=1e-5):
        '''
        return True if are the same (except comments maybe) otherwise False
        '''
        d_k = old['dihedral_coeffs']['coeff_k'] - new['dihedral_coeffs']['coeff_k']
        s_k = old['dihedral_coeffs']['coeff_k'] + new['dihedral_coeffs']['coeff_k']
        if d_k == 0 and s_k == 0:
            ...
        elif abs(d_k / s_k) > float_compare_accuracy:
            return False
        if old['dihedral_coeffs']['coeff_d'] != new['dihedral_coeffs']['coeff_d']:
            return False
        if old['dihedral_coeffs']['coeff_n'] != new['dihedral_coeffs']['coeff_n']:
            return False
        old_label = old['dihedral_coeffs']['comment']
        new_label = new['dihedral_coeffs']['comment']
        if ((old_label and new_label) and (old_label != new_label)):
            return False
        return True

    def util_compare_improper_types(old, new, float_compare_accuracy=1e-5):
        '''
        return True if are the same (except comments maybe) otherwise False
        '''
        d_k = old['improper_coeffs']['coeff_k'] - new['improper_coeffs']['coeff_k']
        s_k = old['improper_coeffs']['coeff_k'] + new['improper_coeffs']['coeff_k']
        if d_k == 0 and s_k == 0:
            ...
        elif abs(d_k / s_k) > float_compare_accuracy:
            return False
        if old['improper_coeffs']['coeff_d'] != new['improper_coeffs']['coeff_d']:
            return False
        if old['improper_coeffs']['coeff_n'] != new['improper_coeffs']['coeff_n']:
            return False
        old_label = old['improper_coeffs']['comment']
        new_label = new['improper_coeffs']['comment']
        if ((old_label and new_label) and (old_label != new_label)):
            return False
        return True

    # merge parameterizations
    mmt_atoms = dfc_mmt.atoms
    mmt_atom_type_ids = set([atom['atom_type_id'] for atom in mmt_atoms])
    mmt_atom_types = [{
       'atom_type_id': atom_type_id,
       'pair_coeffs': None,
       'mass': None,
    } for atom_type_id in mmt_atom_type_ids]
    for mass in dfc_mmt.masses:
        for v in mmt_atom_types:
            if v['atom_type_id'] == mass['mass_id']:
                v['mass'] = mass['mass']
    for pair_coeff in dfc_mmt.pair_coeffs:
        for v in mmt_atom_types:
            if v['atom_type_id'] == pair_coeff['pair_coeff_id']:
                if v['pair_coeffs'] is None:
                    v['pair_coeffs'] = pair_coeff
                else:
                    print('Error: matching pair coeffs')
#    pprint(mmt_atom_types)  # +
    mmt_bonds = dfc_mmt.bonds
    mmt_bond_type_ids = set([bond['bond_type_id'] for bond in mmt_bonds])
    mmt_bond_types = [{
       'bond_type_id': bond_type_id,
       'bond_coeffs': None
    } for bond_type_id in mmt_bond_type_ids]
    for bond_coeff in dfc_mmt.bond_coeffs:
        for v in mmt_bond_types:
            if v['bond_type_id'] == bond_coeff['bond_coeff_id']:
                if v['bond_coeffs'] is None:
                    v['bond_coeffs'] = bond_coeff
                else:
                    print('Error: matching bond coeffs')
#    pprint(mmt_bond_types)  # +
    mmt_angles = dfc_mmt.angles
    if mmt_angles is None:
        mmt_angle_types = []
    else:
        mmt_angle_type_ids = set([angle['angle_type_id'] for angle in mmt_angles])
        mmt_angle_types = [{
           'angle_type_id': angle_type_id,
           'angle_coeffs': None
        } for angle_type_id in mmt_angle_type_ids]
        for angle_coeff in dfc_mmt.angle_coeffs:
            for v in mmt_angle_types:
                if v['angle_type_id'] == angle_coeff['angle_coeff_id']:
                    if v['angle_coeffs'] is None:
                        v['angle_coeffs'] = angle_coeff
                    else:
                        print('Error: matching angle coeffs')
#    pprint(mmt_angle_types)  # +
    mmt_dihedrals = dfc_mmt.dihedrals
    if mmt_dihedrals is None:
        mmt_dihedral_types = []
    else:
        mmt_dihedral_type_ids = set([dihedral['dihedral_type_id']
            for dihedral in mmt_dihedrals])
        mmt_dihedral_types = [{
           'dihedral_type_id': dihedral_type_id,
           'dihedral_coeffs': None
        } for dihedral_type_id in mmt_dihedral_type_ids]
        for dihedral_coeff in dfc_mmt.dihedral_coeffs:
            for v in mmt_dihedral_types:
                if v['dihedral_type_id'] == dihedral_coeff['dihedral_coeff_id']:
                    if v['dihedral_coeffs'] is None:
                        v['dihedral_coeffs'] = dihedral_coeff
                    else:
                        print('Error: matching dihedral coeffs')
#    pprint(mmt_dihedral_types)  # +
    mmt_impropers = dfc_mmt.impropers
    if mmt_impropers is None:
        mmt_improper_types = []
    else:
        mmt_improper_type_ids = set([improper['improper_type_id']
            for improper in mmt_impropers])
        mmt_improper_types = [{
           'improper_type_id': improper_type_id,
           'improper_coeffs': None
        } for improper_type_id in mmt_improper_type_ids]
        for improper_coeff in dfc_mmt.improper_coeffs:
            for v in mmt_improper_types:
                if v['improper_type_id'] == improper_coeff['improper_coeff_id']:
                    if v['improper_coeffs'] is None:
                        v['improper_coeffs'] = improper_coeff
                    else:
                        print('Error: matching improper coeffs')
#    pprint(mmt_improper_types)  # +

    poly_atoms = dfc_poly.atoms
    poly_atom_type_ids = set([atom['atom_type_id'] for atom in poly_atoms])
    poly_atom_types = [{
       'atom_type_id': atom_type_id,
       'pair_coeffs': None,
       'mass': None,
    } for atom_type_id in poly_atom_type_ids]
    for mass in dfc_poly.masses:
        for v in poly_atom_types:
            if v['atom_type_id'] == mass['mass_id']:
                v['mass'] = mass['mass']
    for pair_coeff in dfc_poly.pair_coeffs:
        for v in poly_atom_types:
            if v['atom_type_id'] == pair_coeff['pair_coeff_id']:
                if v['pair_coeffs'] is None:
                    v['pair_coeffs'] = pair_coeff
                else:
                    print('Error: matching pair coeffs')
#    pprint(poly_atom_types)  # +
    poly_bonds = dfc_poly.bonds
    poly_bond_type_ids = set([bond['bond_type_id'] for bond in poly_bonds])
    poly_bond_types = [{
       'bond_type_id': bond_type_id,
       'bond_coeffs': None
    } for bond_type_id in poly_bond_type_ids]
    for bond_coeff in dfc_poly.bond_coeffs:
        for v in poly_bond_types:
            if v['bond_type_id'] == bond_coeff['bond_coeff_id']:
                if v['bond_coeffs'] is None:
                    v['bond_coeffs'] = bond_coeff
                else:
                    print('Error: matching bond coeffs')
#    pprint(poly_bond_types)  # +
    poly_angles = dfc_poly.angles
    if poly_angles is None:
        poly_angle_types = []
    else:
        poly_angle_type_ids = set([angle['angle_type_id']
            for angle in poly_angles])
        poly_angle_types = [{
           'angle_type_id': angle_type_id,
           'angle_coeffs': None
        } for angle_type_id in poly_angle_type_ids]
        for angle_coeff in dfc_poly.angle_coeffs:
            for v in poly_angle_types:
                if v['angle_type_id'] == angle_coeff['angle_coeff_id']:
                    if v['angle_coeffs'] is None:
                        v['angle_coeffs'] = angle_coeff
                    else:
                        print('Error: matching angle coeffs')
#    pprint(poly_angle_types)  # +
    poly_dihedrals = dfc_poly.dihedrals
    if poly_dihedrals is None:
        poly_dihedral_types = []
    else:
        poly_dihedral_type_ids = set([dihedral['dihedral_type_id']
            for dihedral in poly_dihedrals])
        poly_dihedral_types = [{
           'dihedral_type_id': dihedral_type_id,
           'dihedral_coeffs': None
        } for dihedral_type_id in poly_dihedral_type_ids]
        for dihedral_coeff in dfc_poly.dihedral_coeffs:
            for v in poly_dihedral_types:
                if v['dihedral_type_id'] == dihedral_coeff['dihedral_coeff_id']:
                    if v['dihedral_coeffs'] is None:
                        v['dihedral_coeffs'] = dihedral_coeff
                    else:
                        print('Error: matching dihedral coeffs')
#    pprint(poly_dihedral_types)  # +
    poly_impropers = dfc_poly.impropers
    if poly_impropers is None:
        poly_improper_types = []
    else:
        poly_improper_type_ids = set([improper['improper_type_id']
            for improper in poly_impropers])
        poly_improper_types = [{
           'improper_type_id': improper_type_id,
           'improper_coeffs': None
        } for improper_type_id in poly_improper_type_ids]
        for improper_coeff in dfc_poly.improper_coeffs:
            for v in poly_improper_types:
                if v['improper_type_id'] == improper_coeff['improper_coeff_id']:
                    if v['improper_coeffs'] is None:
                        v['improper_coeffs'] = improper_coeff
                    else:
                        print('Error: matching improper coeffs')
#    pprint(poly_improper_types)  # +
    merged_atom_types = copy.deepcopy(mmt_atom_types)
    for atom_type in poly_atom_types:
        is_present = False
        for old_atom_type in merged_atom_types:
            if util_compare_atom_types(old_atom_type, atom_type):
                is_present = True
                break
        if not is_present:
            candidates = set(range(1, len(merged_atom_types) + 2))
            for old_atom_type in merged_atom_types:
                candidates -= set([old_atom_type['atom_type_id']])
            new_atom_type_id = min(candidates)
            merged_atom_types.append({
                'atom_type_id': new_atom_type_id,
                'mass': atom_type['mass'],
                'pair_coeffs': {
                    'coeff_eps': atom_type['pair_coeffs']['coeff_eps'],
                    'coeff_sig': atom_type['pair_coeffs']['coeff_sig'],
                    'comment': atom_type['pair_coeffs']['comment'],
                    'pair_coeff_id': new_atom_type_id}})
    #print('*** mmt ***')
    #for atom_type in mmt_atom_types:
    #    print(atom_type)
    #print('*** poly ***')
    #for atom_type in poly_atom_types:
    #    print(atom_type)
    #print('*** merged ***')
    #for atom_type in merged_atom_types:
    #    print(atom_type)
    # +++++
    merged_bond_types = copy.deepcopy(mmt_bond_types)
    for bond_type in poly_bond_types:
        is_present = False
        for old_bond_type in merged_bond_types:
            if util_compare_bond_types(old_bond_type, bond_type):
                is_present = True
                break
        if not is_present:
            candidates = set(range(1, len(merged_bond_types) + 2))
            for old_bond_type in merged_bond_types:
                candidates -= set([old_bond_type['bond_type_id']])
            new_bond_type_id = min(candidates)
            merged_bond_types.append({
                'bond_type_id': new_bond_type_id,
                'bond_coeffs': {
                    'coeff_k': bond_type['bond_coeffs']['coeff_k'],
                    'coeff_l': bond_type['bond_coeffs']['coeff_l'],
                    'comment': bond_type['bond_coeffs']['comment'],
                    'bond_coeff_id': new_bond_type_id}})
    #print('*** mmt ***')
    #for bond_type in mmt_bond_types:
    #    print(bond_type)
    #print('*** poly ***')
    #for bond_type in poly_bond_types:
    #    print(bond_type)
    #print('*** merged ***')
    #for bond_type in merged_bond_types:
    #    print(bond_type)
    # +++++

    merged_angle_types = copy.deepcopy(mmt_angle_types)
    for angle_type in poly_angle_types:
        is_present = False
        for old_angle_type in merged_angle_types:
            if util_compare_angle_types(old_angle_type, angle_type):
                is_present = True
                break
        if not is_present:
            candidates = set(range(1, len(merged_angle_types) + 2))
            for old_angle_type in merged_angle_types:
                candidates -= set([old_angle_type['angle_type_id']])
            new_angle_type_id = min(candidates)
            merged_angle_types.append({
                'angle_type_id': new_angle_type_id,
                'angle_coeffs': {
                    'coeff_k': angle_type['angle_coeffs']['coeff_k'],
                    'coeff_theta': angle_type['angle_coeffs']['coeff_theta'],
                    'comment': angle_type['angle_coeffs']['comment'],
                    'angle_coeff_id': new_angle_type_id}})
    #print('*** mmt ***')
    #for angle_type in mmt_angle_types:
    #    print(angle_type)
    #print('*** poly ***')
    #for angle_type in poly_angle_types:
    #    print(angle_type)
    #print('*** merged ***')
    #for angle_type in merged_angle_types:
    #    print(angle_type)
    # +++++

    merged_dihedral_types = copy.deepcopy(mmt_dihedral_types)
    for dihedral_type in poly_dihedral_types:
        is_present = False
        for old_dihedral_type in merged_dihedral_types:
            if util_compare_dihedral_types(old_dihedral_type, dihedral_type):
                is_present = True
                break
        if not is_present:
            candidates = set(range(1, len(merged_dihedral_types) + 2))
            for old_dihedral_type in merged_dihedral_types:
                candidates -= set([old_dihedral_type['dihedral_type_id']])
            new_dihedral_type_id = min(candidates)
            merged_dihedral_types.append({
                'dihedral_type_id': new_dihedral_type_id,
                'dihedral_coeffs': {
                    'coeff_k': dihedral_type['dihedral_coeffs']['coeff_k'],
                    'coeff_d': dihedral_type['dihedral_coeffs']['coeff_d'],
                    'coeff_n': dihedral_type['dihedral_coeffs']['coeff_n'],
                    'comment': dihedral_type['dihedral_coeffs']['comment'],
                    'dihedral_coeff_id': new_dihedral_type_id}})
    #print('*** mmt ***')
    #for dihedral_type in mmt_dihedral_types:
    #    print(dihedral_type)
    #print('*** poly ***')
    #for dihedral_type in poly_dihedral_types:
    #    print(dihedral_type)
    #print('*** merged ***')
    #for dihedral_type in merged_dihedral_types:
    #    print(dihedral_type)
    # +++++

    merged_improper_types = copy.deepcopy(mmt_improper_types)
    for improper_type in poly_improper_types:
        is_present = False
        for old_improper_type in merged_improper_types:
            if util_compare_improper_types(old_improper_type, improper_type):
                is_present = True
                break
        if not is_present:
            candidates = set(range(1, len(merged_improper_types) + 2))
            for old_improper_type in merged_improper_types:
                candidates -= set([old_improper_type['improper_type_id']])
            new_improper_type_id = min(candidates)
            merged_improper_types.append({
                'improper_type_id': new_improper_type_id,
                'improper_coeffs': {
                    'coeff_k': improper_type['improper_coeffs']['coeff_k'],
                    'coeff_d': improper_type['improper_coeffs']['coeff_d'],
                    'coeff_n': improper_type['improper_coeffs']['coeff_n'],
                    'comment': improper_type['improper_coeffs']['comment'],
                    'improper_coeff_id': new_improper_type_id}})
    #print('*** mmt ***')
    #for improper_type in mmt_improper_types:
    #    print(improper_type)
    #print('*** poly ***')
    #for improper_type in poly_improper_types:
    #    print(improper_type)
    #print('*** merged ***')
    #for improper_type in merged_improper_types:
    #    print(improper_type)
    # +++++

    dfc_merged = DatafileContent(None)
    dfc_merged.comment = 'merged mmt and polymer'
    dfc_merged.atoms = dfc_mmt.atoms
    dfc_merged.atoms_count = dfc_mmt.atoms_count
    dfc_merged.atom_types = len(merged_atom_types)
    dfc_merged.velocities = None
    dfc_merged.bonds = dfc_mmt.bonds
    dfc_merged.bonds_count = dfc_mmt.bonds_count
    if merged_bond_types is not None:
        dfc_merged.bond_types = len(merged_bond_types)
    else:
        dfc_merged.bond_types = None
    dfc_merged.angles = dfc_mmt.angles
    dfc_merged.angles_count = dfc_mmt.angles_count
    if merged_angle_types is not None:
        dfc_merged.angle_types = len(merged_angle_types)
    else:
        dfc_merged.angle_types = None
    dfc_merged.dihedrals = dfc_mmt.dihedrals
    dfc_merged.dihedrals_count = dfc_mmt.dihedrals_count
    if merged_dihedral_types is not None:
        dfc_merged.dihedral_types = len(merged_dihedral_types)
    else:
        dfc_merged.dihedral_types = None
    dfc_merged.impropers = dfc_mmt.impropers
    dfc_merged.impropers_count = dfc_mmt.impropers_count
    if merged_improper_types is not None:
        dfc_merged.improper_types = len(merged_improper_types)
    else:
        dfc_merged.improper_types = None
    dfc_merged.xlo = dfc_mmt.xlo; dfc_merged.xhi = dfc_mmt.xhi
    dfc_merged.ylo = dfc_mmt.ylo; dfc_merged.yhi = dfc_mmt.yhi
    dfc_merged.zlo = dfc_mmt.zlo; dfc_merged.zhi = dfc_mmt.zhi
    dfc_merged.xy = dfc_mmt.xy
    dfc_merged.xz = dfc_mmt.xz
    dfc_merged.yz = dfc_mmt.yz

    dfc_merged.masses = [{
        'mass_id': atom_type['atom_type_id'],
        'mass': atom_type['mass'],
        'comment': atom_type['pair_coeffs']['comment']
        } for atom_type in merged_atom_types]
    dfc_merged.pair_coeffs = [{
        'pair_coeff_id': atom_type['atom_type_id'],
        'coeff_eps': atom_type['pair_coeffs']['coeff_eps'],
        'coeff_sig': atom_type['pair_coeffs']['coeff_sig'],
        'comment': atom_type['pair_coeffs']['comment']
        } for atom_type in merged_atom_types]
    dfc_merged.bond_coeffs = [{
        'bond_coeff_id': bond_type['bond_type_id'],
        'coeff_k': bond_type['bond_coeffs']['coeff_k'],
        'coeff_l': bond_type['bond_coeffs']['coeff_l'],
        'comment': bond_type['bond_coeffs']['comment']
    } for bond_type in merged_bond_types]
    dfc_merged.angle_coeffs = [{
        'angle_coeff_id': angle_type['angle_type_id'],
        'coeff_k': angle_type['angle_coeffs']['coeff_k'],
        'coeff_theta': angle_type['angle_coeffs']['coeff_theta'],
        'comment': angle_type['angle_coeffs']['comment']
    } for angle_type in merged_angle_types]
    dfc_merged.dihedral_coeffs = [{
        'dihedral_coeff_id': dihedral_type['dihedral_type_id'],
        'coeff_k': dihedral_type['dihedral_coeffs']['coeff_k'],
        'coeff_d': dihedral_type['dihedral_coeffs']['coeff_d'],
        'coeff_n': dihedral_type['dihedral_coeffs']['coeff_n'],
        'comment': dihedral_type['dihedral_coeffs']['comment']
    } for dihedral_type in merged_dihedral_types]
    dfc_merged.improper_coeffs = [{
        'improper_coeff_id': improper_type['improper_type_id'],
        'coeff_k': improper_type['improper_coeffs']['coeff_k'],
        'coeff_d': improper_type['improper_coeffs']['coeff_d'],
        'coeff_n': improper_type['improper_coeffs']['coeff_n'],
        'comment': improper_type['improper_coeffs']['comment']
    } for improper_type in merged_improper_types]

    atom_types_remap = None
    bond_types_remap = None
    angle_types_remap = None
    dihedral_types_remap = None
    improper_types_remap = None

    if mmt_atom_types is not None and len(mmt_atom_types) > 0:
        atom_types_remap = {}
        for poly_atom_type in poly_atom_types:
            for old_atom_type in merged_atom_types:
                if util_compare_atom_types(old_atom_type, poly_atom_type):
                    k = poly_atom_type['atom_type_id']
                    v = old_atom_type['atom_type_id']
                    atom_types_remap[k] = v
    #print(atom_types_remap)  # +
    if mmt_bond_types is not None and len(mmt_bond_types) > 0:
        bond_types_remap = {}
        for poly_bond_type in poly_bond_types:
            for old_bond_type in merged_bond_types:
                if util_compare_bond_types(old_bond_type, poly_bond_type):
                    k = poly_bond_type['bond_type_id']
                    v = old_bond_type['bond_type_id']
                    bond_types_remap[k] = v
    #print(bond_types_remap)  # +
    if mmt_angle_types is not None and len(mmt_angle_types) > 0:
        #print(poly_angle_types)
        angle_types_remap = {}
        for poly_angle_type in poly_angle_types:
            for old_angle_type in merged_angle_types:
                if util_compare_angle_types(old_angle_type, poly_angle_type):
                    k = poly_angle_type['angle_type_id']
                    v = old_angle_type['angle_type_id']
                    angle_types_remap[k] = v
    #print(angle_types_remap)  # +?
    if mmt_dihedral_types is not None and len(mmt_dihedral_types) > 0:
        dihedral_types_remap = {}
        for poly_dihedral_type in poly_dihedral_types:
            for old_dihedral_type in merged_dihedral_types:
                if util_compare_dihedral_types(old_dihedral_type,
                    poly_dihedral_type):
                    k = poly_dihedral_type['dihedral_type_id']
                    v = old_dihedral_type['dihedral_type_id']
                    dihedral_types_remap[k] = v
    #print(dihedral_types_remap)  # +?
    if mmt_improper_types is not None and len(mmt_improper_types) > 0:
        improper_types_remap = {}
        for poly_improper_type in poly_improper_types:
            for old_improper_type in merged_improper_types:
                if util_compare_improper_types(old_improper_type,
                    poly_improper_type):
                    k = poly_improper_type['improper_type_id']
                    v = old_improper_type['improper_type_id']
                    improper_types_remap[k] = v
    #print(improper_types_remap)  # +?

    #print(atom_types_remap)  # +?
    #print(bond_types_remap)  # +?
    #print(angle_types_remap)  # +?
    #print(dihedral_types_remap)  # +?
    #print(improper_types_remap)  # +?

    # merge geometries
    bbox_mmt = get_bbox(dfc_mmt)
    bbox_poly = get_bbox(dfc_poly)
    mmt_lx = dfc_mmt.xhi - dfc_mmt.xlo
    mmt_ly = dfc_mmt.yhi - dfc_mmt.ylo
    mmt_lz = dfc_mmt.zhi - dfc_mmt.zlo
    poly_cell_lx = dfc_poly.xhi - dfc_poly.xlo
    poly_cell_ly = dfc_poly.yhi - dfc_poly.ylo
    poly_cell_lz = dfc_poly.zhi - dfc_poly.zlo
    poly_bbox_lx = bbox_poly['bbox_xhi'] - bbox_poly['bbox_xlo']
    poly_bbox_ly = bbox_poly['bbox_yhi'] - bbox_poly['bbox_ylo']
    poly_bbox_lz = bbox_poly['bbox_zhi'] - bbox_poly['bbox_zlo']
    dxlo = bbox_mmt['bbox_xlo'] - bbox_poly['bbox_xlo']
    dylo = bbox_mmt['bbox_ylo'] - bbox_poly['bbox_ylo']
    dzlo = bbox_mmt['bbox_zlo'] - bbox_poly['bbox_zlo']
    dxhi = bbox_mmt['bbox_xhi'] - bbox_poly['bbox_xhi']
    dyhi = bbox_mmt['bbox_yhi'] - bbox_poly['bbox_yhi']
    dzhi = bbox_mmt['bbox_zhi'] - bbox_poly['bbox_zhi']
    mmt_cell_lo_bbox_lo = dfc_mmt.zlo - bbox_mmt['bbox_zlo']
    mmt_cell_hi_bbox_hi = bbox_mmt['bbox_zhi'] - dfc_mmt.zhi
    bbox_mmt_empty = {
        'bbox_xlo': bbox_mmt['bbox_xlo'], 'bbox_xhi': bbox_mmt['bbox_xhi'],
        'bbox_ylo': bbox_mmt['bbox_ylo'], 'bbox_yhi': bbox_mmt['bbox_yhi'],
    }
    if mmt_cell_lo_bbox_lo > mmt_cell_hi_bbox_hi:
        bbox_mmt_empty['bbox_zlo'] = bbox_mmt['bbox_zhi']
        bbox_mmt_empty['bbox_zhi'] = dfc_mmt.zhi
    else:
        bbox_mmt_empty['bbox_zlo'] = dfc_mmt.zlo
        bbox_mmt_empty['bbox_zhi'] = bbox_mmt['bbox_zlo']
    if mmt_lx < poly_bbox_lx:
        print('error! small box in x direction!')
        print('mmt bbox lx:', mmt_lx)
        print('poly lx:', poly_bbox_lx)
        return False
    if mmt_ly < poly_bbox_ly:
        print('error! small box in y direction!')
        print('mmt ly', mmt_ly)
        print('poly bbox ly:', poly_bbox_ly)
        return False
    if mmt_lz < poly_bbox_lz:
        print('error! small box in z direction!')
        print('mmt lz:', mmt_lz)
        print('poly bbox lz:', poly_bbox_lz)
        return False
    empty_lz = bbox_mmt_empty['bbox_zhi'] - bbox_mmt_empty['bbox_zlo']
    # how much polys can be inserted along axes
    #print(empty_lz, poly_bbox_lz)
    n_poly_x = int(mmt_lx // (poly_bbox_lx + 2*border))
    n_poly_y = int(mmt_ly // (poly_bbox_ly + 2*border))
    n_poly_z = int(empty_lz // (poly_bbox_lz + 2*border))
    print(n_poly_x, n_poly_y, n_poly_z)
    print('so, inserting {0} polymers'.format(n_poly_x * n_poly_y * n_poly_z))
    for idx_x in range(n_poly_x):
        tmp_xlo = (bbox_mmt_empty['bbox_xlo']
            + (2*border + poly_bbox_lx) * idx_x)
        for idx_y in range(n_poly_y):
            tmp_ylo = (bbox_mmt_empty['bbox_ylo']
                + (2*border + poly_bbox_ly) * idx_y)
            for idx_z in range(n_poly_z):
                tmp_zlo = (bbox_mmt_empty['bbox_zlo']
                    + (2*border + poly_bbox_lz) * idx_z)
                # insert single poly molecule now
                tmp_atom_idx = len(dfc_merged.atoms) + 1
                try:
                    tmp_bond_idx = len(dfc_merged.bonds) + 1
                except TypeError:
                    tmp_bond_idx = 1
                try:
                    tmp_angle_idx = len(dfc_merged.angles) + 1
                except TypeError:
                    tmp_angle_idx = 1
                try:
                    tmp_dihedral_idx = len(dfc_merged.dihedrals) + 1
                except TypeError:
                    tmp_dihedral_idx = 1
                try:
                    tmp_improper_idx = len(dfc_merged.impropers) + 1
                except TypeError:
                    tmp_improper_idx = 1
                tmp_molecule_tag = max(atom['molecule-tag']
                    for atom in dfc_merged.atoms) + 1

                atoms_delta = len(dfc_merged.atoms)
                # TODO
                for atom in dfc_poly.atoms:
                    try:
                        new_atom_type = atom_types_remap[atom['atom_type_id']]
                    except TypeError:
                        new_atom_type = atom['atom_type_id']
                    dfc_merged.atoms.append({
#                        'atom_id': tmp_atom_idx,
                        'atom_id': atoms_delta + atom['atom_id'],
                        'molecule-tag': tmp_molecule_tag,
                        'atom_type_id': new_atom_type,
                        'q': atom['q'],
                        'x': (bbox_mmt_empty['bbox_xlo']
                              + atom['x'] - bbox_poly['bbox_xlo']
                              + border + idx_x * (poly_bbox_lz + 2*border)),
                        'y': (bbox_mmt_empty['bbox_ylo']
                              + atom['y'] - bbox_poly['bbox_ylo']
                              + border + idx_y * (poly_bbox_lz + 2*border)),
                        'z': (bbox_mmt_empty['bbox_zlo']
                              + atom['z'] - bbox_poly['bbox_zlo']
                              + border + idx_z * (poly_bbox_lz + 2*border)),
                        'nx': atom['nx'],
                        'ny': atom['ny'],
                        'nz': atom['nz'],
                        'comment': atom['comment']
                    })
                    new_z = (bbox_mmt_empty['bbox_zlo']
                             + atom['z'] - bbox_poly['bbox_zlo']
                             + border + idx_z * (poly_cell_lz + 2*border))
                    #assert(dfc_merged.zlo < new_z < dfc_merged.zhi)

                    tmp_atom_idx += 1
                dfc_merged.atoms_count = len(dfc_merged.atoms)

                if dfc_merged.bonds is None:
                    dfc_merged.bonds = []
                for bond in dfc_poly.bonds:
                    #print(bond)
                    try:
                        new_bond_type = bond_types_remap[bond['bond_type_id']]
                    except TypeError:
                        new_bond_type = bond['bond_type_id']
                    dfc_merged.bonds.append({
                        'bond_id': tmp_bond_idx,
                        'bond_type_id': new_bond_type,
                        'atom_one_id': bond['atom_one_id'] + atoms_delta,
                        'atom_two_id': bond['atom_two_id'] + atoms_delta
                    })
                    tmp_bond_idx += 1
                dfc_merged.bonds_count = len(dfc_merged.bonds)

                if dfc_merged.angles is None:
                    dfc_merged.angles = []
                for angle in dfc_poly.angles:
                    #print(angle)
                    try:
                        new_angle_type = angle_types_remap[angle['angle_type_id']]
                    except TypeError:
                        new_angle_type = angle['angle_type_id']
                    dfc_merged.angles.append({
                        'angle_id': tmp_angle_idx,
                        'angle_type_id': new_angle_type,
                        'atom_one_id': angle['atom_one_id'] + atoms_delta,
                        'atom_two_id': angle['atom_two_id'] + atoms_delta,
                        'atom_three_id': angle['atom_three_id'] + atoms_delta
                    })
                    tmp_angle_idx += 1
                dfc_merged.angles_count = len(dfc_merged.angles)

                if dfc_merged.dihedrals is None:
                    dfc_merged.dihedrals = []
                for dihedral in dfc_poly.dihedrals:
                    #print(dihedral)
                    try:
                        new_dihedral_type = \
                            dihedral_types_remap[dihedral['dihedral_type_id']]
                    except TypeError:
                        new_dihedral_type = dihedral['dihedral_type_id']
                    dfc_merged.dihedrals.append({
                        'dihedral_id': tmp_dihedral_idx,
                        'dihedral_type_id': new_dihedral_type,
                        'atom_one_id': dihedral['atom_one_id'] + atoms_delta,
                        'atom_two_id': dihedral['atom_two_id'] + atoms_delta,
                        'atom_three_id': dihedral['atom_three_id'] + atoms_delta,
                        'atom_four_id': dihedral['atom_four_id'] + atoms_delta
                    })
                    tmp_dihedral_idx += 1
                dfc_merged.dihedrals_count = len(dfc_merged.dihedrals)

                if dfc_merged.impropers is None:
                    dfc_merged.impropers = []
                for improper in dfc_poly.impropers:
                    #print(improper)
                    try:
                        new_improper_type = \
                            improper_types_remap[improper['improper_type_id']]
                    except TypeError:
                        new_improper_type = improper['improper_type_id']
                    dfc_merged.impropers.append({
                        'improper_id': tmp_improper_idx,
                        'improper_type_id': new_improper_type,
                        'atom_one_id': improper['atom_one_id'] + atoms_delta,
                        'atom_two_id': improper['atom_two_id'] + atoms_delta,
                        'atom_three_id': improper['atom_three_id'] + atoms_delta,
                        'atom_four_id': improper['atom_four_id'] + atoms_delta
                    })
                    tmp_improper_idx += 1
                dfc_merged.impropers_count = len(dfc_merged.impropers)
    return dfc_merged


if __name__ == '__main__':
    dezired_nz = 10

    fname_mmt = 'data_structures/mmt.data'
    fname_mmt_mod = 'data_structures/mmt_and_mod.data'
    fname_poly = 'data_structures/pa6_x10.data'
    dfc_mmt = DatafileContent(fname_mmt)
    dfc_mmt_mod = DatafileContent(fname_mmt_mod)
    dfc_poly = DatafileContent(fname_poly)

    bbox_poly = get_bbox(dfc_poly)

    dfc_mmt_mod.zhi += dezired_nz * (bbox_poly['bbox_zhi'] - bbox_poly['bbox_zlo'])

    dfc_merged = insert_poly_into_mmt(dfc_mmt_mod, dfc_poly)

    dfc_merged.write('new_mmt.data')
