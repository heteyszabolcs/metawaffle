__author__ = 'sgalan'
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def check(input_file, outdir):
    def simple_plot(ax):
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

    peaks_dict = defaultdict(list)
    widths = []
    with open(input_file, 'r') as r:
        for line in r:
            c, s, e = line.split()
            widths.append(abs(int(s)-int(e)))
            peaks_dict[c].append(int((int(s) + int(e)) / 2))

    distances = []
    for c in peaks_dict:
        ordered = sorted(peaks_dict[c])
        for i in range(len(ordered)-1):
            distances.append(np.log10(ordered[i+1] - ordered[i]))


    plt.ioff()
    plt.switch_backend('Agg')
    fig, ax = plt.subplots(2, 1)
    ax[0].set_title('Peak evaluation')
    _ = ax[0].hist(distances,bins=100,histtype='stepfilled',color='Grey')
    ax[0].set_xlabel('Distance contiguous peaks, log$^{10}$(bp)')
    simple_plot(ax[0])
    _ = ax[1].hist(widths,bins=25, histtype='stepfilled',color='black')
    ax[1].set_xlabel('Width (bp)')
    simple_plot(ax[1])
    plt.tight_layout()
    plt.savefig(outdir+'test_evaluation.png', dpi=100)
    plt.close(fig)
