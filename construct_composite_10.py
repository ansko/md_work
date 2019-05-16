import random

from datafile_content import DatafileContent
from datafile_content_modifiers import get_bbox
from insert_poly_into_mmt import insert_poly_into_mmt

from utils import move_some_atoms


if __name__ == '__main__':
    '''
    Construct composite based on mix.data
        3480 atoms, like in old segregated system
    '''

    # mmt and mod extruded from mix.data
    fname_mmt_mod = 'data_structures/mmt_and_mod.data'
    # some relaxed poly molecule (not from mix.data!)
    fname_poly = 'data_structures/pa6_x10.data'

    dezired_nz = 50

    dfc_mmt_mod = DatafileContent(fname_mmt_mod)
    dfc_poly = DatafileContent(fname_poly)

    # correct charges in mmt: obts -> obos
    q_obts = -1.1688  # net charge =  0.5424
    q_obos = -1.1808  # net charge = -0.0336 (much less!)
    for atom in  dfc_mmt_mod.atoms:
        if atom['comment'] == 'obts':
            atom['comment'] = 'obos';
            atom['q'] = q_obos

    # extend box:
    bbox_poly = get_bbox(dfc_poly)
    dz = dezired_nz * (bbox_poly['bbox_zhi'] - bbox_poly['bbox_zlo'])
    dfc_mmt_mod.zhi += dz

    # insert polymers
    dfc_merged = insert_poly_into_mmt(dfc_mmt_mod, dfc_poly)

    # move half of modifiers
    dfc_mmt_mod.zhi += 40  # additional space for modifiers
    chosen_modifiers_idcs = set()
    while len(chosen_modifiers_idcs) != 6:
        chosen_modifiers_idcs.add(random.choice(range(12)))
    chosen_modifiers_idcs = list(sorted(chosen_modifiers_idcs))
    for modifier_idx in chosen_modifiers_idcs:
        atoms_idcs = range(721 + modifier_idx*70, 791 + modifier_idx*70)
        move_some_atoms(dfc_mmt_mod, atoms_idcs, dz=dz+40)
    dfc_merged.zhi += 40

    # save data
    dfc_merged.write('new_mmt.data')
