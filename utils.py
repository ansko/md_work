def move_some_atoms(dfc, atoms_idcs_list, dx=0, dy=0, dz=0):
    for atom in dfc.atoms:
        if atom['atom_id'] in atoms_idcs_list:
            atom['x'] += dx
            atom['y'] += dy
            atom['z'] += dz


def change_atoms(dfc, condition, remap):
    done = 0
    for atom in dfc.atoms:
        try:
            atom['type'] = remap[condition(atom)]
            done += 1
        except KeyError:
            pass
    return done


if __name__ == '__main__':
   from datafile_content import DatafileContent

   dfc = DatafileContent('na_mmt.125000.data')

   # into atom labels:
   condition = lambda atom: str(atom['q'])
   remap = {
       '1.575':      'Al',
       '1.36':       'Mg',
       '2.1':        'Si',
       '-1.05':       'O',
       '-1.0808':     'O',
       '-0.95':       'O',
       '-1.1688':     'O',
       '0.4245':      'H',
       '1.0':         'Na'
   }
   print(change_atoms(dfc, condition, remap))

   dfc.write('000.data')
