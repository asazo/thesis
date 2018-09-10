#!/usr/bin/env python3


import sys
import os
sys.path.append("/home/asazo/repos/grains-2017-2018/code")
import pandas as pd
from ggstats.statistics import hist_areas
import numpy as np

# Load results from imagej
folder = sys.argv[1]
listdir = os.listdir(folder)
areas = []
for file in listdir:
    areas.append(pd.read_csv(folder+"/"+file)['Area'].as_matrix())
areas = np.hstack(areas)
print("Loaded %d areas." % len(areas))

hist_areas("Test", areas, "areas", ".", use_data=False, use_mean=True)