#usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import statistics

#load in the data
data = pd.read_csv("../data/SNV_free_energies.tsv", sep = '\t')

#compute the stats
free_energies = np.array(data["Free_Energy_(kcal/mol)"])
mean = np.mean(free_energies)
median = np.median(free_energies)
free_energy_wt = 338.445
standard_deviation = statistics.stdev(free_energies)

#calculate ninety fifth confidence interval
ci= stats.norm.interval(confidence=0.95, loc=mean, scale=stats.sem(free_energies))
print(ci)

#write the stats to a txt file
with open('../outputs/free_energy_stats.txt', 'w') as f:
    f.write("free energy wt:  338.445 kcal_mol")
    f.write("mean:  " + str(mean))
    f.write("median:  " + str(median))
    f.write("standard deviation:  " + str(standard_deviation))
    f.write("95th confidence interval:  " + str(ci))


#make the histogram
hist, bins = np.histogram(free_energies)
plt.figure()
logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]))
plt.hist(data, bins=logbins, ec = 'k')
plt.plot(mean)
plt.plot(median)
plt.plot(free_energy_wt)
plt.xlabel("Free Energy")
plt.ylabel("Frequency")
plt.title("Frequency of Free Energy Values for SNV Missense Mutations")
plt.show()
