#!/usr/bin/env python3

import sys
import os
sys.path.append("/home/asazo/repos/grains-2017-2018/code")
import pandas as pd
from ggstats.statistics import hist_areas
import ggstats
ggstats.config.set_latex(True)
ggstats.config.set_fancy_output(True)
import numpy as np

# Load results from imagej, they are stored in each grain*/cuts folder
#listdir = [l for l in os.listdir('.') if os.path.isdir(os.path.join('.',l)) and l[:6] == "grains"]
listdir = ['grains1', 'grains2', 'grains3']
areas = []
for folder in listdir:
    cuts_folder = os.path.join(folder, "cuts")
    print("Reading folder %s" %  cuts_folder)
    filename_x = os.path.join(cuts_folder, "results-x.xls")
    filename_y = os.path.join(cuts_folder, "results-y.xls")
    filename_z = os.path.join(cuts_folder, "results-z.xls")
    try:
        areas_x = pd.read_csv(filename_x, sep='\t')['Area'].values
        areas.append(areas_x)
    except:
        print("Couldn't load %s" % filename_x)
    try:
        areas_y = pd.read_csv(filename_y, sep='\t')['Area'].values
        areas.append(areas_y)
    except:
        print("Couldn't load %s" % filename_y)
    try:
        areas_z = pd.read_csv(filename_z, sep='\t')['Area'].values
        areas.append(areas_z)
    except:
        print("Couldn't load %s" % filename_z)

    #areas.append(pd.read_csv(folder+"/"+file)['Area'].as_matrix())
areas = np.sort(np.hstack(areas))
areas = areas[areas > 4]
print("Max grain area is %.16f" % np.max(areas))
print("Loaded %d areas." % len(areas))
hist_areas("", areas, "areas", ".", use_data=True, use_mean=True, bins=50)