#usr/bin/env python3

#import packages
import numpy as np
import pandas as pd

#load in the data
somatic = pd.read_csv("../data/tissue_hotspots_somatic.tsv", sep = '\t')
free_energies = pd.read_csv("../data/SNV_free_energies.tsv", sep = '\t')
ins_free_energies = pd.read_csv("../data/SNV_ins_free_energies.tsv", sep = '\t')
del_free_energies = pd.read_csv("../data/SNV_del_free_energies.tsv", sep = '\t')

#make the big lists
possible_muts_column = []
possible_free_energies_column = []
possible_msa_column = []

#big loop
for codon in somatic["Codon"]:
    possible_muts = []
    possible_free_energies = []
    possible_msas = []
    #SNV missense
    for mutation in free_energies["Mutation"]:
        if (str(codon) in str(mutation)) == True:
            possible_muts.append(mutation)
            mini_df = free_energies.loc[free_energies["Mutation"] == mutation]
            for free_energy in mini_df["Free_Energy_(kcal/mol)"]:
                possible_free_energies.append(free_energy)
            for msa in mini_df["MSA_Name"]:
                possible_msas.append(msa)
        else:
            pass

    #SNV insertion
    for insertion in ins_free_energies["Mutation"]:
        if (str(codon) in str(insertion)) == True:
            possible_muts.append(insertion)
            mini_ins_df = ins_free_energies.loc[ins_free_energies["Mutation"] == insertion]
            for ins_free_energy in mini_ins_df["Free_Energy_(kcal/mol)"]:
                possible_free_energies.append(ins_free_energy)
            for ins_msa in mini_ins_df["MSA_Name"]:
                possible_msas.append(ins_msa)
        else:
            pass

    #SNV deletion
    for deletion in del_free_energies["Mutation"]:
        if (str(codon) in str(deletion)) == True:
            possible_muts.append(deletion)
            mini_del_df = del_free_energies.loc[del_free_energies["Mutation"] == deletion]
            for del_free_energy in mini_del_df["Free_Energy_(kcal/mol)"]:
                possible_free_energies.append(del_free_energy)
            for del_msa in mini_del_df["MSA_Name"]:
                possible_msas.append(del_msa)
        else:
            pass
    
    possible_muts_column.append(possible_muts)
    possible_free_energies_column.append(possible_free_energies)
    possible_msa_column.append(possible_msas)

#add columns to dataframe
somatic["Possible_Mutations"] = possible_muts_column
somatic["Free_Energy_(kcal/mol)"] = possible_free_energies_column
somatic["MSA_Name"] = possible_msa_column
somatic.to_csv("../outputs/updated_somatic.tsv", sep = '\t')
