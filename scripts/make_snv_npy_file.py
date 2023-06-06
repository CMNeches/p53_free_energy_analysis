#usr/bin/env python3

#import packages
import numpy as np
import pandas as pd

#load in the data
data = pd.read_csv('../data/SNV_free_energies.tsv', sep = '\t')["Mutation"]
array = np.array(data)

#make the npy file
np.save("../outputs/SNV_free_energy.npy", array)
