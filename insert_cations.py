import random
import time


from datafile_content import DatafileContent


def insert_cations(dfc_mmt, N, **ff):
    ff_type = ff['ff_type']
    ff_mass = ff['ff_mass']
    ff_charge = ff['ff_charge']
    ff_pair_coeffs = ff['ff_pair_coeffs']

    done = 0
    xlo = dfc_mmt.xlo; xhi = dfc_mmt.xhi; lx = xhi - xlo
    ylo = dfc_mmt.ylo; yhi = dfc_mmt.yhi; ly = yhi - ylo
    zlo = dfc_mmt.zlo; zhi = dfc_mmt.zhi; lz = zhi - zlo
    random.seed(int(time.time()))
    too_close_distance_sq = 5
    new_atom_type = max([atom['atom_type_id'] for atom in dfc_mmt.atoms]) + 1
    while done < N:
        coeff_x = random.random()
        coeff_y = random.random()
        coeff_z = random.random()
        x = xlo + lx * coeff_x
        y = ylo + ly * coeff_y
        z = zlo + lz * coeff_z
        is_close = False
        for atom in dfc_mmt.atoms:
            dx = atom['x'] - x
            dy = atom['y'] - y
            dz = atom['z'] - z
            if dx**2 + dy**2 + dz**2 < too_close_distance_sq:
                is_close = True
                break
        if is_close:
            continue
        else:
            dfc_mmt.atoms.append({
                'atom_id': len(dfc_mmt.atoms) + 1,
                'molecule-tag': 1,
                'atom_type_id': new_atom_type,
                'q': ff_charge,
                'x': x,
                'y': y,
                'z': z,
                'nx': 0,
                'ny': 0,
                'nz': 0,
                'comment': ff_type
            })
        done += 1
    dfc_mmt.atom_types += 1
    dfc_mmt.atoms_count += N
    dfc_mmt.masses.append({
        'mass_id': new_atom_type,
        'mass': ff_mass,
        'comment': ff_type
    })
    dfc_mmt.pair_coeffs.append({
        'pair_coeff_id': new_atom_type,
        'coeff_eps': ff_pair_coeffs['eps'],
        'coeff_sig': ff_pair_coeffs['sig'],
        'comment': ff_type
    })

if __name__ == '__main__':
    dfc_mmt = DatafileContent('data_structures/mmt.data')
    q = sum([atom['q'] for atom in dfc_mmt.atoms])
    ff_na = {
        'ff_type': 'Na',
        'ff_mass': '22.99',
        'ff_charge': 1.0,
        'ff_pair_coeffs': {
            'eps': 0.13009998714677298,
            'sig': 2.3500126639309213
        }
    }
    insert_cations(dfc_mmt, 12, **ff_na)
    dfc_mmt.write('na_mmt.data')
