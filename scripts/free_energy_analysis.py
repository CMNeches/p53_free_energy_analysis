#usr/bin/env python3

#import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import regex

#load in the data
germline_hotspots = pd.read_csv("../data/tissue_hotspots_germline.tsv", sep = '\t')
somatic_hotspots = pd.read_csv("../data/tissue_hotspots_somatic.tsv", sep = '\t')
free_energies = pd.read_csv("../data/SNV_free_energies.tsv", sep = '\t')

#set up the hotspots
germline_hotspots = germline_hotspots["Codon"].to_list()
somatic_hotspots = somatic_hotspots["Codon"].to_list()

#separate the free energies
codons = free_energies["Mutation"].str.replace('\D','', regex = True).to_list()
final_codons = []
for codon in codons:
    final_codon = int(codon)
    final_codons.append(final_codon)
free_energies = free_energies["Free_Energy_(kcal/mol)"].to_list()

#determine the free energies of hotspots:
germline_free_energies = []
somatic_free_energies = []

#determine if the codon in free energies is a hotspot, if so, pull the free energy
#instead of only pulling first instance, pull each instance and take the means of their free_energies
germline_hotspots_to_graph = []
somatic_hotspots_to_graph = []
def free_energy_puller(mutation_type):
    for codon in final_codons:
        if codon in mutation_type:
            count = int(mutation_type.count(codon))
            index = final_codons.index(codon)
            free_energy = free_energies[index]
            if mutation_type == germline_hotspots:
                germline_free_energies.append(free_energy)
                germline_hotspots_to_graph.append(codon)
            else:
                somatic_free_energies.append(free_energy)
                somatic_hotspots_to_graph.append(codon)
    return(germline_free_energies)
    return(somatic_free_energies)
    return(germline_hotspots_to_graph)
    return(somatic_hotspots_to_graph)

free_energy_puller(germline_hotspots)
free_energy_puller(somatic_hotspots)

#save the files
germline_df = pd.DataFrame(germline_free_energies)
somatic_df = pd.DataFrame(somatic_free_energies)

germline_df.to_csv("../outputs/germline_free_energies.tsv", sep = '\t')
somatic_df.to_csv("../outputs/somatic_free_energies.tsv", sep = '\t')

updated_germline_hotspots = pd.DataFrame(germline_hotspots_to_graph)
updated_somatic_hotspots = pd.DataFrame(somatic_hotspots_to_graph)

updated_germline_hotspots.to_csv("../outputs/germline_hotspots_to_graph.tsv", sep = '\t')
updated_somatic_hotspots.to_csv("../outputs/somatic_hotspots_to_graph.tsv", sep = '\t')
