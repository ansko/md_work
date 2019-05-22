'''
MMT
'''


def mmt_density_profile(dfc, multip):
    zmin = zmax = dfc.atoms[0]['z']
    density_profile = {}
    for atom in dfc.atoms:
        k = int(multip * atom['z'])
        try:
            density_profile[k] += 1
        except KeyError:
            density_profile[k] = 1
        zmin = min(zmin, atom['z'])
        zmax = max(zmax, atom['z'])
    lz = zmax - zmin
    for atom in dfc.atoms:
        k = int(multip * (atom['z'] + lz))
        try:
            density_profile[k] += 1
        except KeyError:
            density_profile[k] = 1
    for atom in dfc.atoms:
        k = int(multip * (atom['z'] - lz))
        try:
            density_profile[k] += 1
        except KeyError:
            density_profile[k] = 1

    print('lz:', zmax - zmin, zmin, zmax)

    return density_profile


if __name__ == '__main__':
    fname = 'na_mmt.230000.data'
    result = mmt_density_profile([fname])
