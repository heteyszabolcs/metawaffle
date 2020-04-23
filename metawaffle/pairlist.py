from collections import OrderedDict, defaultdict
from itertools   import combinations
from os          import path
import math
import logging
logger = logging.getLogger('')

def binning_bed(peak_file, windows_span, max_dist, outdir,
                name, chrom_sizes, windows, **kwargs):
    '''Input bed file of ChIP peaks and bin into desired resolution of Hi-C'''

    def read_line_feature(line):
        '''
        Get information per peak of a feature +/-
        '''
        c, p1, p2, f = line.split()[:4]
        return c, int(round((int(p1) + int(p2)) / float(2))) , f

    def read_line_no_feature(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, int(round((int(p1) + int(p2)) / float(2))) , ''

    peaks = open(peak_file, "r")

    # findout if bed file contain features, or only coordinates
    line = peaks.next()
    try:
        read_line_feature(line)
        read_line = read_line_feature
    except ValueError:
        read_line = read_line_no_feature
    peaks.seek(0)

    bin_coordinate = set((c, p, f) for c, p, f in map(read_line, peaks)
                          if p > windows_span)  # take into account to add windows span both 
    logger.info('[MAIN]: Doing all the combinations...')
    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features
    
    wsp = (windows_span * 2)
    mdr = max_dist

    pairs = ((a, b) for a, b in combinations(bin_coordinate, 2) if a[0] == b[0]
             and wsp <= abs(b[1] - a[1]) <= mdr)

    logger.info('[MAIN]: Writing pairlists.')

    intervals = {}
    for (chromosome1, bs1, f1), (chromosome2,  bs2, f2) in pairs:
        distance = abs(bs1 - bs2)
        if bs1 < bs2:
            for lower, upper in windows:
                if lower < distance <= upper:
                    intervals.setdefault((lower, upper), {})
                    intervals[(lower, upper)].setdefault(
                        (f1, f2), []).append((chromosome1, bs1, chromosome2, bs2))
        if bs1 > bs2:
            for lower, upper in windows:
                if lower < distance <= upper:
                    intervals.setdefault((lower, upper), {})
                    intervals[(lower, upper)].setdefault((f2,f1), []).append((chromosome2, bs2, chromosome1, bs1))
    
    # define categories and write pairs
    for beg, end in intervals:
        logger.info('[MAIN]: Writing interval: %i, %i'%(beg, end))
        for f1, f2 in intervals[(beg, end)]:
            w = open(path.join(outdir, '%s_%d_%d.tsv' % (
                name, beg, end)), 'w')
            for c1, s1, c2, s2 in intervals[(beg, end)][(f1, f2)]:
                # check chromosome length
                
                new_start1, new_end1 = s1 - windows_span, s1 + windows_span
                new_start2, new_end2 = s2 - windows_span, s2 + windows_span
                    
                if new_end1 <= chrom_sizes[c1] and new_end2 <= chrom_sizes[c1]:
                    w.write('%s:%d-%d\t%s:%d-%d\n' % (
                        c1, new_start1, new_end1, c2, new_start2, new_end2))
            w.close()
