#This script is for copying the convergence logs from the nice run on Hamilton, currently ongoing as of the 23rd July.
#
#Should be easy enough -- I only have half an hour until the train...


import numpy as np
import matplotlib.pyplot as plt
import os, sys
import matplotlib

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 12,        # Default font size
    "axes.labelsize": 12,
    "axes.titlesize": 8,
    "xtick.labelsize": 12,
    "ytick.labelsize": 8,
})

fig_width = 443.57848/72

cs = matplotlib.colormaps.get_cmap('tab10')

if False:
    #Copy files across
    for i in range(8):
        string = 'scp vgjn10@hamilton8.dur.ac.uk:/nobackup/vgjn10/projects/vsw_forecasting/data/cma_data/%d_log.csv ./data/cma_data/' % i
        os.system(string)

#Establish plot baseline and things
def make_nicetitle(id):
    if (id//2)%2 == 0:
        is_pfss = True
        model = "PFSS"
    else:
        is_pfss = False
        model = "Outflow"
    if (id%2) == 0:
        rss = 2.5
    else:
        rss = 5.0
    if (id//4) == 0:
        source = "GONG"
    else:
        source = "HMI"

    rss_string = "r_{ss}"
    nicetitle = f"{model}, ${rss_string} = {rss}$, {source}"

    return nicetitle
#Target persistence score:
target_score = 0.6731977470239495
original_scores = [1.39500977, 1.33923572, 1.22588864, 1.07992013, 1.07768663, 1.30469469, 1.05301729, 0.99063013]

fig = plt.figure(figsize = (fig_width,0.6*fig_width))
batch_cadence = 16
for i in range(8):
    fname = './paper/data/cma_data/%d_log.csv' %i
    data = np.loadtxt(fname, delimiter = ',')
    #Get convergence score
    scores = data[:,1]
    #Remove outliers
    scores[scores > 5] = np.nan

    batch_scores = []; batch_centres = []
    nbatches = int(len(scores)/batch_cadence)
    for n in range(nbatches):
        if not np.all(np.isnan(scores[batch_cadence*n:batch_cadence*(n+1)])):
            batch_scores.append(np.nanmin(scores[batch_cadence*n:batch_cadence*(n+1)]))
            batch_centres.append(0.5*(batch_cadence*n + batch_cadence*(n+1)))
    plt.plot(batch_centres, batch_scores, c = cs(i), linewidth=1.0, label = make_nicetitle(i))
    plt.plot(scores, c = cs(i), linewidth=0.1, alpha=0.1)

    plt.axhline(original_scores[i], c = cs(i), linestyle='dashed', linewidth=1.0)

plt.axhline(target_score, c = 'black', linestyle='dashed', label='Persistence')
handles, labels = plt.gca().get_legend_handles_labels()
fig.legend(handles, labels,
    loc="lower center",
    ncol=3,                  # adjust as needed
    bbox_to_anchor=(0.5, -0.0), fontsize=8)

plt.yscale('log')
plt.xlim(-10,150)
plt.ylim(0.6,3.0)
plt.gca().set_yticks([], minor=True)
plt.gca().set_yticks([0.75,1.0,1.5,2.5])
plt.gca().set_yticklabels([0.75,1.0,1.5,2.5])
plt.ylabel('Combined Skill Score')
plt.xlabel('Iteration')
plt.tight_layout(rect=[0, 0.13, 1, 1])
plt.savefig('./paper/plots/6_plot_cma_convergence.pdf')
plt.show()



