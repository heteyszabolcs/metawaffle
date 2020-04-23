from collections import OrderedDict
import subprocess
import uuid
import csv
import glob
import os
from collections import defaultdict
from cPickle import dump,load
import random
import numpy as np
import logging
logger = logging.getLogger('')

def extract_coordinates(peak_list,resolution,section_pos,tmpdir,chromosome,badcols, names):
    '''Chunk file into multiple, and write them in parallel per file write coord of 10,000 peak pairs
    Write a dictionary depending of pairs of peak per target'''

    unique_filename = str(uuid.uuid4())
    w = open(tmpdir+chromosome+'_{}.tmp'.format(unique_filename),'wa')
    position = 0
    for line in peak_list:
        peak1, peak2 = line.split()
        chr1, beg1, end1 = peak1.replace(':','-').split('-')
        chr2, beg2, end2 = peak2.replace(':','-').split('-')
        pos = section_pos[chr1][0]
        label = chr1+'_'+str(position)
        names[label] = (peak1, peak2)
        if int(beg1) < int(beg2):
            start_bin1 = pos + (int(beg1) / resolution)
            end_bin1 = pos + (int(end1) / resolution)

            start_bin2 = pos + (int(beg2) / resolution)
            end_bin2 = pos + (int(end2) / resolution)
        else:
            start_bin1 = pos + (int(beg2) / resolution)
            end_bin1 = pos + (int(end2) / resolution)

            start_bin2 = pos + (int(beg1) / resolution)
            end_bin2 = pos + (int(end1) / resolution)
        for x, p1 in enumerate(xrange(start_bin1, end_bin1+1)):
            for y, p2 in enumerate(xrange(start_bin2, end_bin2+1)):
                if p1 in badcols or p2 in badcols:
                    continue
                w.write('{}\t{}\t{}\t{}\t{}\n'.format(p1, p2, x, y, label))
        position += 1
    w.close()

def random_line(fname):
    '''Extract random line from pairs peak file'''
    lines = open(fname).read().splitlines()
    return random.choice(lines)

def eq_pos(pos1, pos2):
    return pos1 == pos2

def greater_pos(pos1, pos2):
    return pos1 > pos2

def readfiles(file1,file2,chromosome, avg_nrm):
    def split_line1(l):
        if len(l.split()) == 4:
            a, b, c, d = l.split()
            return (int(a), int(b)), d
        else:
            a, b, c = l.split()
            return (int(a), int(b)), c

    def split_line2(l):
        a, b, c, d, f = l.split()
        return (int(a), int(b)), int(c), int(d), f
    logger.info('[MAIN]: Reading matrix and peaks...')
    fh1 = open(file1)
    fh2 = open(file2)
    pos1, nrm = split_line1(fh1.next())
    pos2, x, y, f = split_line2(fh2.next())
    try:
        while True:
            if eq_pos(pos1, pos2):
                avg_nrm[f][(x,y)] = float(nrm)
                pos2_ = pos2
                pos2, x, y, f = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated but different cell in metamatrix
                    pos1, nrm = split_line1(fh1.next())
            elif greater_pos(pos1, pos2):
                pos2, x, y, f = split_line2(fh2.next())
            else:
                pos1, nrm = split_line1(fh1.next())
                
    except StopIteration:
        fh1.close()
        fh2.close()
    logger.info('[MAIN]: Finished.')

def write_matrices(avg_nrm, outdir, name, size, names):
    '''To obtain the mean matrix, divide raw and norm per passages'''
    logger.info('[MAIN]: Writing submatrices...')
    w = open(outdir+'/matrices_%s.tsv'%(name),'wa')
    total_array = [[] for i in range(len(avg_nrm.keys()))]
    for n, label in enumerate(avg_nrm.keys()):
        array_nrm = np.zeros((size, size))
        for coordinates in avg_nrm[label]:
            array_nrm[coordinates] = avg_nrm[label][coordinates]
        # concatenate the rows in 1D
        flat_list = np.ndarray.tolist(array_nrm.flatten())
        total_array[n] = flat_list
        total_array[n].insert(0,names[label][0])
        total_array[n].insert(1, names[label][1])
    # write results to txt file (row-wise)
    wr = csv.writer(w)
    wr.writerows(total_array)
    w.close()
