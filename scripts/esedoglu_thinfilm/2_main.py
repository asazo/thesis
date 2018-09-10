#!/usr/bin/env python3

"""
    Generate automatic cuts over a given 3D grain structure
"""

import numpy as np
import sys
import os
import scipy.io
import pandas as pd
from skimage import measure, morphology
import pymesh
from sklearn.feature_extraction.image import img_to_graph
import pickle
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.cm as cm

dims = (128,128,128)


def load_grain(grains, k):
    """
        Load single grain from matlab source.
        The grain is a surface in 3D.
    """
    grain = -np.ones(dims)
    ind = grains[k][0]-1
    [x, y, z] = np.unravel_index(ind, dims, order='F')
    val = grains[k][1]
    grain[y,x,z] = val
    verts, faces = measure.marching_cubes_classic(grain, 0, spacing=(1, 1, 1))
    return verts, faces


def build_plane(direction, level, limit, eps=2):
    """
        Generate the plane to perform the cut
        given a certain direction x, y or z
        and a certain level between limits
    """
    if direction == 'y':
        Xa = np.array([limit, level, limit])
        Xb = np.array([0, level+eps, 0])
    elif direction == 'z':
        Xa = np.array([limit, limit, level])
        Xb = np.array([0, 0, level+eps])
    elif direction == 'x':
        Xa = np.array([level, limit, limit])
        Xb = np.array([level+eps, 0, 0])
    return pymesh.generate_box_mesh(Xa, Xb)


def generate_cuts(grains, n_grains, direction='x', levels=np.array([64]), zoom=4):
    """
        Generate the cuts that can be plotted to image to be analyzed.
    """
    # Column to be rescued, depends on direction
    if direction == 'x':
        cols = [0,1,2]
    elif direction == 'y':
        cols = [1,0,2]
    elif direction == 'z':
        cols = [2,0,1]
    results_per_level = [[None for i in range(n_grains)] for level in levels]
    # For each grain intersect with each of the different planes
    for i in range(n_grains):
        # Load grain from MATLAB output
        #g = load_grain(i)
        # Get vertices and faces (grain triangulation)
        #verts, faces = build_surface(g)
        verts, faces = grains[i]
        # Build high level mesh
        mesh = pymesh.form_mesh(verts, faces)
        # Iterate over planes
        for j, level in enumerate(levels):
            plane = build_plane(direction, level, 127)
            result = pymesh.boolean(mesh, plane, operation="intersection", engine="carve")
            if len(result.vertices) > 0:
                results_per_level[j][i] = result
    # Build cuts
    all_cuts = [None for results in results_per_level]
    for i, (level, results) in enumerate(zip(levels, results_per_level)):
        finals = []
        # Filter results, some Nones are just non intersecting grains
        # The result is only PyMesh.Mesh objects
        results = [r for r in results if r is not None]
        #  a single cut over a level is the union of each PyMesh.Mesh objects
        for r in results:
            mask = r.vertices[:,cols[0]] == level
            tmp = r.vertices[mask]
            finals.append(tmp[:, [cols[1],cols[2]]])
        finals = np.vstack(finals)
        # Amplify image by a factor given by zoom
        finals = (finals * zoom).astype(int)
        # Binarize output
        M = np.zeros((np.max(finals)+1, np.max(finals)+1))
        M[finals[:,1], finals[:,0]] = 1
        # Dilate and skeletonize
        dit = morphology.binary_dilation(M, selem=np.ones((4,4)))
        skt = morphology.skeletonize(dit)
        all_cuts[i] = skt
    return all_cuts


def save_cuts(all_cuts, grains_folder, direction='x'):
    """
        Generate a single image
    """
    # Border size in pixels
    border_size = 20
    # Number of cuts
    n_cuts = len(all_cuts)
    # Original shape of cuts
    shape = all_cuts[0].shape
    shape = np.array([shape[1], shape[0]])
    new_shape = shape
    # Shape of output image
    new_shape[1] = shape[1] + border_size
    new_shape[1] = new_shape[1] * n_cuts
    # Generate empty image with shape (shape + border_size) * n_cuts
    img = np.zeros(new_shape)
    print("Output image shape",new_shape)
    for i, cut in enumerate(all_cuts):
        #print(i*(new_shape[0]+border_size), (i+1)*shape[0] + i*border_size)
        img[:, i*(new_shape[0]+border_size):(i+1)*shape[0] + i*border_size] = cut[::-1]
        #plt.imshow(cut[::-1], cmap=cm.binary)
    fig = plt.figure(figsize=(n_cuts*10, 10))
    plt.imshow(img, cmap=cm.binary)
    plt.axis("off")
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
    plt.savefig(grains_folder+"/%s.png" % (direction), box_inches='tight', dpi=100,)
    plt.close()
    plt.clf()


def help():
    print("./2_main.py FILE")

def main(argv):
    # Load the output of gbm3d.m after many iterations.
    if len(argv) != 2:
        help()
        return
    matlab_source = argv[1]
    grains = scipy.io.loadmat(matlab_source)['grains']
    n_grains = len(grains)
    print("Loaded %d grains" % n_grains)
    # Create folder structure
    base_folder = os.path.join("./", matlab_source.split("/")[-1][:-4])
    grains_folder = base_folder + "/grains"
    cuts_folder = base_folder + "/cuts"

    # Check if the files were saved before
    if os.path.isdir(base_folder):
        filenames = os.listdir(grains_folder)
        filenames = [grains_folder+"/"+filename for filename in filenames]
        filenames.sort()
        print("Loading saved grains...")
        grains = [None for i in filenames]
        for i, filename in enumerate(tqdm(filenames)):
            with open(filename, 'rb') as f:
                ret = pickle.load(f)
                grains[i] = ret
        print("Loaded.")
        levels = np.arange(10,120,20, dtype=int)
        print("Generating %d levels" % len(levels))
        for dimension in ['x', 'y', 'z']:
            all_cuts = generate_cuts(grains, n_grains, dimension, levels)
            save_cuts(all_cuts, cuts_folder, dimension)
    else:
        os.mkdir(base_folder)
        os.mkdir(grains_folder)
        os.mkdir(cuts_folder)
        print("Generating grains and exporting to %s..." % grains_folder)
        for g_id in tqdm(range(n_grains)):
            ret = load_grain(grains, g_id)
            with open(grains_folder + "/%d.pkl" % g_id, 'wb') as f:
                pickle.dump(ret, f)
        print("Done.")


if __name__ == '__main__':
    main(sys.argv)