#Delete this script once I'm done with it

import os

plot_fnames = os.listdir('./plots/overfit_converge')


print(plot_fnames)

for fname in plot_fnames:
    if fname[:3] == 'bad':
        print(fname[4:], len(fname))
        number = fname[4:-4]
        new_fname = './plots/overfit_converge/bad_%03d.png' % int(number)
        print(fname, new_fname)
        string = f'mv ./plots/overfit_converge/{fname} {new_fname}'
        print(string)
        os.system(f'mv ./plots/overfit_converge/{fname} {new_fname}')
