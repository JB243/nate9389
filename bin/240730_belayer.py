import os
from pprint import pprint
from csv import reader
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, homogeneity_score,v_measure_score
from sklearn import linear_model,preprocessing
from scipy.optimize import minimize, Bounds
from scipy.linalg import orth, eigh
from scipy.stats import mode
import statsmodels.api
import statsmodels as sm
import anndata
import scanpy as sc
from glmpca import glmpca
import csv
from tqdm import trange
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from importlib import reload

import sys
sys.path.append('src') # Load Belayer functions
from utils_IO import *
from harmonic import *
from region_cost_fun import *
from dprelated import *
from dp_post_processing import *
from general_helper_funcs import * # Cong added the helper function script from belayer-develop to belayer repo
from precompute_likelihood import *
from dp_linear_boundary import *
import svg
from slideseq_helpers import *

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Example usage: python app.py --path data')
    parser.add_argument('--path', type=str, help='Input path (e.g., ../240523_adata_NC.h5ad)')
    args = parser.parse_args()

    adata = sc.read_h5ad(args.path)
    gene_labels=np.array(adata.var.index)
    count=adata.X.copy().toarray().T
    G,N=count.shape
    coords=adata.obs[['array_row', 'array_col']].values
    glmpca_res = glmpca.glmpca(count, 2*7, fam="poi", penalty=10, verbose=True)
    F_glmpca = glmpca_res['factors']
    
    max_nlayers=8 # max number of layers
    rotation_angle_list=[i for i in range(360)]

    loss_array,label_dict=rotation_dp(F_glmpca.T, coords, Lmax=max_nlayers, use_buckets=True, 
                          num_buckets=150,rotation_angle_list=rotation_angle_list)

    min_loss = float('inf')
    min_loss_index = -1
    
    for idx, layer_losses in enumerate(loss_array):
        min_layer_loss = min(layer_losses)
        if min_layer_loss < min_loss:
            min_loss = min_layer_loss
            min_loss_index = idx

    best_angle=rotation_angle_list[min_loss_index]
    best_L=2
    belayer_labels, belayer_rotated_coords=visualized_rotated_belayer_output(label_dict, coords, best_L, 
                                                                         best_angle,
                                                   rotation_angle_list=rotation_angle_list)
    adata.obs['belayer'] = belayer_labels
    sc.pl.spatial(adata, color = 'belayer', save = args.path.split('/')[-1].split('.')[0] + '.png')
    adata.write_h5ad('belayered_' + args.path.split('/')[-1])
