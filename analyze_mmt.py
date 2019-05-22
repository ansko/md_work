import matplotlib
import matplotlib.pyplot as plt

from structural_tools.datafile_content import DatafileContent
from analyzers.mmt_density_profile import mmt_density_profile


if __name__ == '__main__':
    multip = 10

    fnames = ['na_mmt.{0}000.data'.format(i) for i in range(275, 380, 5)]

    kss = set()
    profiles = []
    for fname in fnames:
        profiles.append(mmt_density_profile(DatafileContent(fname), multip))

    values = {}
    for profile in profiles:
        for k, v in profile.items():
            try:
                values[k] += v
            except KeyError:
                values[k] = v

    x_label_str = r'z'
    y_label_str = r'Atoms count'
    out_fname = '1_average_modifier_distance.pdf'

    fig = plt.figure()
    plt.xlabel(x_label_str)
    plt.ylabel(y_label_str)

    result = {}
    for k in values.keys():
        try:
            result[k] += values[k]
        except KeyError:
            result[k] = values[k]

    plt.plot([k / multip - 15.5 for k in sorted(result.keys())],
        [result[k] for k in sorted(result.keys())], color='k')
    plt.gca().axes.get_yaxis().set_visible(False)

    #plt.gca().set_xlim([15.5, 25.5])
    plt.gca().set_xlim([0, 10])

    fig.savefig(out_fname)

    print('done!')
