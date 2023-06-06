#usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics

#import the data
germline_free_energies = pd.read_csv("../outputs/germline_free_energies.tsv", sep = '\t')
germline_free_energies = germline_free_energies["0"].to_list()
somatic_free_energies = pd.read_csv("../outputs/somatic_free_energies.tsv", sep = '\t')
somatic_free_energies = somatic_free_energies["0"].to_list()
germline_hotspots_for_histogram = pd.read_csv("../outputs/germline_hotspots_to_graph.tsv", sep = '\t')
germline_hotspots_for_histogram = germline_hotspots_for_histogram["0"].to_list()
somatic_hotspots_for_histogram = pd.read_csv("../outputs/somatic_hotspots_to_graph.tsv", sep = '\t')
somatic_hotspots_for_histogram = somatic_hotspots_for_histogram["0"].to_list()

#make the histograms

#germline
plt.figure()
plt.hist(germline_free_energies, ec = 'k')
plt.title("Germline Hotspot Free Energies (kcal/mol)")
plt.xlabel("Frequency")
plt.gca().set_xscale("log")
plt.ylabel("Free Energy (kcal/mol)")
plt.savefig("../outputs/germline_hotspot_free_energy_freqs.png")

plt.figure()
plt.bar(germline_hotspots_for_histogram, germline_free_energies, ec = 'k')
plt.title("Germline Hotspot Codons vs. Free Energy (kcal/mol)")
plt.xlabel("Germline Hotspot Codon")
plt.gca().set_xscale("log")
plt.ylabel("Free Energy (kcal/mol)")
plt.savefig("../outputs/germline_hotspot_codon_free_energy_codon_graph.png")

plt.scatter(germline_hotspots_for_histogram, germline_free_energies)
plt.savefig("../outputs/germline_hotspot_scatter.png")
#somatic
plt.figure()
plt.hist(somatic_free_energies, ec = 'k')
plt.title("Somatic Hotspot Free Energies")
plt.xlabel("Frequency")
plt.gca().set_xscale("log")
plt.ylabel("Free Energy (kcal/mol")
plt.savefig("../outputs/somatic_hotspot_free_energy_freqs.png")

plt.figure()
plt.hist(somatic_hotspots_for_histogram, somatic_free_energies, ec = 'k')
plt.title("Somatic Hotspot Codons vs. Free Energy (kcal/mol)")
plt.xlabel("Somatic Hotspot Codon")
plt.gca().set_xscale("log")
plt.ylabel("Free Energy (kcal/mol)")
plt.savefig("../outputs/somatic_hotspot_codon_free_energy_codon_graph.png")
