#!/usr/bin/env python

from skimage import measure
import mayavi as mlab
import numpy as np
import sys

# Notice that indices here goes from 0 to N-1
dims = (128,128,128)

def load_grain(k):
    u = -np.ones(dims)
    ind = grains[k][0]-1
    [x, y, z] = np.unravel_index(ind, dims, order='F')
    val = grains[k][1]
    u[y,x,z] = val
    return u

def build_surface(grain):
    verts, faces = measure.marching_cubes_classic(grain, 0, spacing=(1, 1, 1))
    return verts, faces

def plot_grain_mayavi(grain):
    verts, faces = build_surface(grain)
    mlab.triangular_mesh([vert[0] for vert in verts],
                         [vert[1] for vert in verts],
                         [vert[2] for vert in verts],
                          faces)
    mlab.show()

if __name__ == '__main__':
    try:
        idx = int(sys.argv[1])
    except:
        print("%s is not an integer" % sys.argv[1])
        bad = True
    if not bad:
        grain = load_grain(idx)
        plot_grain_mayavi(grain)