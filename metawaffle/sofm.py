from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import scipy as sp
import glob
import scipy.ndimage
import matplotlib.gridspec as gridspec
from neupy import architectures, algorithms, environment
from sklearn.preprocessing import scale
import pickle

def matrix_to_line(matrix):
    line = defaultdict(list)
    n = 0
    for x in range(len(matrix)):
        for y in range(len(matrix)):
            line[n] = matrix[x,y]
            n += 1
    return line.values()

def draw_grid(sofm, matrices, output_features, coords, outdir, size):
    data = np.asarray(matrices)
    clusters = sofm.predict(output_features).argmax(axis=1)
    grid_height, grid_weight = sofm.features_grid
    info_array = np.zeros((grid_height, grid_weight))

    for row_id in range(grid_height):
        for col_id in range(grid_weight):
            index = row_id * grid_height + col_id
            w = open(outdir+'Cluster_%s.bed'%(str(index)), 'wa')

            clustered_samples = data[clusters == index]
            coords_cluster = coords[clusters == index]
            if len(clustered_samples) > 0:
                sample = [clustered_samples[0][x:x+size] for x in range(0, len(clustered_samples[0]), size)]
                info_array[row_id, col_id] = len(clustered_samples)
                # write keep coordinates per cluster
                for x in coords_cluster:
                     w.write(x+'\n')
            else:
                # If we don't have samples in cluster then
                # it means that there is a gap in space
                sample = np.zeros((size, size))
            w.close()
    return sample, info_array


def plot_sofm(pairs_file, outdir, size, grid, label_name):
    all_clusters_files = glob.glob(outdir+'Cluster_*.bed')
    pairs_cluster = defaultdict(int)
    cluster_pairs = defaultdict(list)
    for f in all_clusters_files:
        cluster = int(f.split('/')[-1].split('.')[0].split('_')[-1])
        with open(f,'r') as r:
            for line in r:
                pair1 = '-'.join(line.split('\n')[0].split('-')[0:2])
                pair2 ='-'.join(line.split('\n')[0].split('-')[2:])
                pairs_cluster[(pair1,pair2)] = cluster
                cluster_pairs[cluster].append((pair1, pair2))
    total_cells = float(size * size)
    #number_clusters = len(cluster_pairs.keys())
    number_clusters = len(all_clusters_files)
    matrix_info = np.zeros((number_clusters,size*size))
    length_cluster = defaultdict(int)
    pairs_mat = defaultdict(list)
    with open(pairs_file,'r') as r:
        for line in r:
            line = line.strip().split(',')
            pair1 = line[0]
            pair2 = line[1]
            cluster = pairs_cluster[(pair1,pair2)]
            matrix = map(float, line[2:])
            pairs_mat[(pair1, pair2)] = matrix
            m = np.asanyarray(matrix).reshape((size, size))
            matrix_info[cluster] += matrix
            length_cluster[cluster] += 1
    matrix_mean = []
    for c in range(number_clusters):
        m = matrix_info[c]
        matrix = np.asarray(m).reshape((size,size))
        matrix_mean.append(np.divide(matrix, float(length_cluster[c])))
    
    supermatrix= np.zeros((size*grid, size*grid))
    min_num = 2.
    max_index = grid * grid
    idx_dict = defaultdict(list)
    for row_id in range(grid):
        for col_id in range(grid):
            index = row_id * grid + col_id
            if length_cluster[index] > min_num:
                idx_dict[index] = [row_id,col_id]
    for index in idx_dict:
        row, col = idx_dict[index]
        supermatrix[row*size:(row*size)+size, col*size:(col*size)+size] = matrix_mean[index]
    plt.ioff()
    plt.switch_backend('Agg')
    fig = plt.figure(figsize=(4,4))
    im = plt.imshow(np.log2(supermatrix), cmap='RdBu_r', interpolation='none', vmax=2, vmin=-2)
    plt.axis('off')

    vals = np.arange(0, len(supermatrix)+size, size)
    _ = [plt.axvline(i-0.5, color='gray', linewidth=0.2) for i in vals]
    _ = [plt.axhline(i-0.5, color='gray', linewidth=0.2) for i in vals]
    plt.tight_layout()
    plt.savefig(outdir+'%s_SOFM.png'%(label_name), dpi=1000)
    plt.close(fig)
