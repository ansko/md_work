import random
import sys

from datafile_content import DatafileContent
from datafile_content_modifiers import get_bbox
from insert_poly_into_mmt import insert_poly_into_mmt
from utils import move_some_atoms


if __name__ == '__main__':
    '''
    Constructs ternary nanocomposite == clay + modifier + polymer
    Polymer is 10-mer or 25-mer or 50-mer or 100-mer
    (polymerization degree == sys.argv[1]).
    Always monomers number in the system == 7500
    (for various polymer molecules count).

    One function instead of four.
    '''

    polymerization_degree = int(sys.argv[1])

    # some relaxed poly molecule (not from mix.data!)
    if polymerization_degree == 10:
        dezired_n = 750
        fname_poly = 'data_structures/pa6_x10.data'
        dezired_nz = dezired_n / 3 / 5
        out_fname = 'comp_mon10_n750.data'
    elif polymerization_degree == 25:
        dezired_n = 300
        fname_poly = 'data_structures/pa6_x25.data'
        dezired_nz = dezired_n / 1 / 2
        out_fname = 'comp_mon25_n300.data'
    elif polymerization_degree == 50:
        dezired_n = 150
        fname_poly = 'data_structures/pa6_x50.data'
        dezired_nz = dezired_n / 1 / 2
        out_fname = 'comp_mon50_n150.data'
    elif polymerization_degree == 100:
        dezired_n = 75
        fname_poly = 'data_structures/pa6_x100.data'
        dezired_nz = 75
        out_fname = 'comp_mon100_n75.data'
    else:
        print(polymerization_degree, 'is not supported')
        sys.exit()

    # mmt and mod extruded from mix.data
    fname_mmt_mod = 'data_structures/mmt_mod_331.data'
    dfc_mmt_mod = DatafileContent(fname_mmt_mod)

    # correct charges in mmt: obts -> obos
    q_obts = -1.1688  # net charge =  0.5424
    q_obos = -1.1808  # net charge = -0.0336 (much less!)
    for atom in  dfc_mmt_mod.atoms:
        if atom['comment'] == 'obts':
            atom['comment'] = 'obos';
            atom['q'] = q_obos

    dfc_poly = DatafileContent(fname_poly)

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
    dfc_merged.write(out_fname)

    print('done!')
