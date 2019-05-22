import matplotlib
import matplotlib.pyplot as plt

from structural_tools.datafile_content import DatafileContent
from analyzers.mmt_density_profile import mmt_density_profile


if __name__ == '__main__':
    multip = 10

    fnames_300 = ['na_mmt_300_diffusion/na_mmt.{0}000.data'.format(i)
        for i in range(501, 751)]
    fnames_500 = ['na_mmt_500_diffusion/na_mmt.{0}000.data'.format(i)
        for i in range(501, 751)]

    # 721 -- 732
    old_coords = {
        atom['atom_id']: {
            'x': atom['x'],
            'y': atom['y']
        } for atom in DatafileContent(fnames_300[0]).atoms[720:732]
    }

    Ns300 = [0]
    for fname in fnames_300:
        print(fname)

        dfc = DatafileContent(fname)
        new_coords = {
            atom['atom_id']: {
                'x': atom['x'],
                'y': atom['y']
            } for atom in dfc.atoms[720:732]
        }
        lx = dfc.xhi - dfc.xlo
        ly = dfc.yhi - dfc.ylo

        n = Ns300[-1]
        for k in new_coords.keys():
            dx = abs(new_coords[k]['x'] - old_coords[k]['x'])
            dy = abs(new_coords[k]['y'] - old_coords[k]['y'])
            dx = min(dx, abs(lx - dx))
            dy = min(dy, abs(ly - dy))
            if dx > 1 or dy > 1:
                n += 1
                #print(k, dx, dy)

        Ns300.append(n)
        old_coords = new_coords

    ##########

    old_coords = {
        atom['atom_id']: {
            'x': atom['x'],
            'y': atom['y']
        } for atom in DatafileContent(fnames_500[0]).atoms[720:732]
    }

    Ns500 = [0]
    for fname in fnames_500:
        print(fname)

        dfc = DatafileContent(fname)
        new_coords = {
            atom['atom_id']: {
                'x': atom['x'],
                'y': atom['y']
            } for atom in dfc.atoms[720:732]
        }
        lx = dfc.xhi - dfc.xlo
        ly = dfc.yhi - dfc.ylo

        n = Ns500[-1]
        for k in new_coords.keys():
            dx = abs(new_coords[k]['x'] - old_coords[k]['x'])
            dy = abs(new_coords[k]['y'] - old_coords[k]['y'])
            dx = min(dx, abs(lx - dx))
            dy = min(dy, abs(ly - dy))
            if dx > 1 or dy > 1:
                n += 1
                #print(k, dx, dy)

        Ns500.append(n)
        old_coords = new_coords

    #print(Ns)

    fig = plt.figure()
    plt.xlabel('Время, фс')
    plt.ylabel('Количество перескоков')

    line300, = plt.plot(range(len(Ns300)), Ns300, color='k')
    line500, = plt.plot(range(len(Ns500)), Ns500, color='r')

    matplotlib.pyplot.legend([line300, line500], ['300 К', '500 К'],
        loc="upper left")

    fig.savefig('jumps.pdf')

    print('done!')
