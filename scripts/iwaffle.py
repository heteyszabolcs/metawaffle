from argparse                     import ArgumentParser
import pysam
from collections import OrderedDict
from datetime import datetime
import subprocess
import uuid
import csv
import glob
import os
from collections import defaultdict
from cPickle import dump,load
import random
import numpy as np
import math

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
            start_bin1 = pos + int(round((int(beg1) / float(resolution))))
            end_bin1 = pos + int(round((int(end1) / float(resolution))))

            start_bin2 = pos + int(round((int(beg2) / float(resolution))))
            end_bin2 = pos + int(round((int(end2) / float(resolution))))
        else:
            start_bin1 = pos + int(round((int(beg2) / float(resolution))))
            end_bin1 = pos + int(round((int(end2) / float(resolution))))

            start_bin2 = pos + int(round((int(beg1) / float(resolution))))
            end_bin2 = pos + int(round((int(end1) / float(resolution))))
        for x, p1 in enumerate(xrange(start_bin1, end_bin1)):
            for y, p2 in enumerate(xrange(start_bin2, end_bin2)):
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

def readfiles(file1,file2,chromosome):
    def split_line1(l):
        a, b, c, d = l.split()
        return (int(a), int(b)), c, d
    def split_line2(l):
        a, b, c, d, f = l.split()
        return (int(a), int(b)), int(c), int(d), f

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Reading BAM and peaks ...'
    fh1 = open(file1)
    fh2 = open(file2)
    pos1, raw, nrm = split_line1(fh1.next())
    pos2, x, y, f = split_line2(fh2.next())
    try:
        while True:
            if eq_pos(pos1, pos2):
                avg_nrm[f][(x,y)] = float(nrm)
                pos2_ = pos2
                pos2, x, y, f = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated but different cell in metamatrix
                    pos1, raw, nrm = split_line1(fh1.next())
            elif greater_pos(pos1, pos2):
                pos2, x, y, f = split_line2(fh2.next())
            else:
                pos1, raw, nrm = split_line1(fh1.next())

    except StopIteration:
        fh1.close()
        fh2.close()
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Finished'

def write_matrices(avg_nrm, outdir, name, size, names):
    '''To obtain the mean matrix, divide raw and norm per passages, plot if wanted'''
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Writing submatrices ...'
    w = open(outdir+'matrices_%s.tsv'%(name),'wa')
    total_array = [[] for i in range(len(avg_nrm.keys()))]
    print size
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

def main():
    opts = get_options()
    inbam        = opts.inbam
    resolution   = opts.resolution
    pairs_file    = opts.pairs_file
    tmpdir       = opts.tmpdir
    outdir       = opts.outdir
    ncpus        = opts.ncpus
    name         = opts.name
    biases       = opts.biases
    mats         = opts.mats
    sample       = opts.sample
    windows_span = opts.windows_span

    if opts.get_all:
        get_all = opts.get_all
    else:
        get_all = False

    bamfile = pysam.AlignmentFile(inbam,'rb')
    sections = OrderedDict(zip(bamfile.references,[x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    if opts.get_all == False:
        print  datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Getting randomly: ',sample, 'from ', name
        #split file  peaks per chromosome

        random_file = open(outdir+name+'_'+str(sample)+'rnd.tsv','wa')
        for n in range(sample):
            line = random_line(pairs_file)
            random_file.write('{}\n'.format(line))
        random_file.close()

        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Splitting peak pairs per chromosome...'
        fh = subprocess.Popen("awk -F: '{print >> " + '"' +  tmpdir + '"' +  " $1; close($1)}' %s " % ((outdir+name+'_'+str(sample)+'rnd.tsv')),
                                  shell=True)
        fh.communicate() # wait until the process finishes

    else:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Getting all submatrices from ', name
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Splitting peak pairs per chromosome...'
        fh = subprocess.Popen("awk -F: '{print >> " + '"' +  tmpdir + '"' +  " $1; close($1)}' %s " % ((pairs_file)),
                                  shell=True)
        fh.communicate() # wait until the process finishes

    chromosomes_file = glob.glob(tmpdir+'*')
    badcols = load(open(biases))['badcol']
    global avg_nrm
    avg_nrm = defaultdict(lambda: defaultdict(float))
    global names
    names = defaultdict(tuple)

    for peak_list in chromosomes_file:
        chromosome = peak_list.split('/')[-1]
        peakfile = open(peak_list,'r')
        # parallel job to write coordinates, splitting peak file, write tmp files
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Writing coordinates files & sorting...', chromosome
        extract_coordinates(peak_list=peakfile, chromosome=chromosome,resolution=resolution,
                                   section_pos=section_pos,tmpdir=tmpdir,badcols=badcols, names=names)

        tmp_chr = glob.glob(tmpdir+'%s_*'%chromosome)
        fh = subprocess.Popen(("sort -k1,2n --parallel=24 -S 20%% %s > %s%s_sorted")%(tmp_chr[0],
                                                                tmpdir,chromosome),shell=True)
        fh.communicate()
        os.system("rm "+tmp_chr[0])

        file1 = mats+'%s_bam_%ikb.tsv'%(chromosome,resolution/1000)
        file2 = '%s%s_sorted'%(tmpdir,chromosome)
        check_size = open(file2,'r')
        nsize = check_size.readlines()
        if len(nsize) > 0:
            readfiles(file1,file2,chromosome)
#            os.system("rm %s%s_sorted"%(tmpdir,chromosome))
#            os.system("rm %s%s"%(tmpdir,chromosome))
        else:
            print 'No information for ', chromosome, ' all badcolumns :('
            os.system("rm %s%s_sorted"%(tmpdir, chromosome))
            os.system("rm %s%s"%(tmpdir, chromosome))
            
    size = ((windows_span * 2) / resolution)
    if len(avg_nrm) != 0:
        write_matrices(avg_nrm, outdir, name, size, names) # write individual matrix to 1D line

    else:
        print 'No information in this subsample: ', name



def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    parser.add_argument('-i','--pairs', dest='pairs_file',required=True, default=False,
                        help='''Pairs peaks file''')
    parser.add_argument('-bam','--bam',dest='inbam',required=True, default=False,
                        help= 'Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-t','--tmp',dest='tmpdir',required=True, default=False,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-o','--outdir',dest='outdir',default=True,help='output directory')
    parser.add_argument('-C','--cpus',dest='ncpus',type=int,default=8,
                        help='''[%(default)s] number of cpus to be used for parsing the HiC-BAM file''')
    parser.add_argument('-n','--name',dest='name',default=True, help = 'Output name')
    parser.add_argument('-b','--biases',dest='biases',default=True, help = 'Biases')
    parser.add_argument('-mats','--mats',dest='mats',default=True, help = 'Folder where matrices are located')
    parser.add_argument('-sample','--sample',dest='sample',default=True,type=int,
                        help = 'Number of pairs that will be extracted randomly', required=False)
    parser.add_argument('-s','--windows_span',dest='windows_span',default=True,type=int,
                        help = 'bp added around target location')
    parser.add_argument('-A','--get_all',dest='get_all',default=False,type=bool, required=False,
                        help = 'To get all the submatrices')

    opts = parser.parse_args()

    return opts

if __name__=='__main__':
    exit(main())
