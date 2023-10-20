#motivo
from Bio import SeqIO
import csv
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import statistics
import random as rd
import re
import itertools
import pandas as pd
from Bio.Seq import Seq
from itertools import product
####################################################################################################################################################################################


#CHECK WHETHER MOTOIF IS IN A SEQUENCE
def generate_non_unique_combinations(motif):
    non_unique_bases = {
        'R': 'AG',
        'Y': 'CT',
        'W': 'AT',
        "N": "ACGT"
    }
    motifs_with_combinations = [motif]
    for code, bases in non_unique_bases.items():
        if code in motif:
            motif_combinations = [motif.replace(code, replacement) for replacement in bases]
            motifs_with_combinations.extend(motif_combinations)
    return motifs_with_combinations

def check_motifs_in_sequence(sequence, motifs, operator="OR"):
    sequence = Seq(sequence)
    occurrence_counter = 0
    for motif in motifs:
        motif_combinations = generate_non_unique_combinations(motif)
        for motif_combination in motif_combinations:
            forward_found = motif_combination in sequence
            reverse_found = motif_combination in sequence.reverse_complement()
            if forward_found or reverse_found:
                occurrence_counter += 1
    if operator == "AND" and occurrence_counter == len(motifs):
        return True
    if operator == "OR" and occurrence_counter > 0:
        return True
    return False

# FIND MOST SIGNIFICANT MOTIFS FOR GENE LIST
def motif_finder(bed_file, fasta_file, gene_ids, csv_file, k=6, size=30):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    significant_motifs = []
    nucleotides = ["A", "C", "G", "T"]
    all_kmers = []
    # for k in range(5, k + 1):
    # all_kmers.extend(kmers_of_length_k)
    all_kmers = [''.join(mer) for mer in itertools.product(nucleotides, repeat=k)]
    print("working with",len(all_kmers),"k-mers")
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end), gene_name))
            gene_list.append(gene_name)
    for binding_site in all_kmers:
        if all_kmers.index(binding_site) % 10 == 0:
            progress_message = "Progress: " + str(all_kmers.index(binding_site) / len(all_kmers) * 100)
            print(progress_message, end="\r")
        myc_targets = []
        search_range = 50
        for chromosome, start, end, gene_name in tss_data:
            sequence = genome[chromosome][int(start) - search_range: int(end) + search_range]
            if check_motifs_in_sequence(sequence, [binding_site], operator="OR"):
                myc_targets.append(gene_name)
        observed_matches = []
        for i in range(len(gene_ids) - size + 1):
            genes_subset = gene_ids[i:i + size]
            matches = len(set(genes_subset) & set(myc_targets))
            observed_matches.append(matches)
        if len(observed_matches) >= len(gene_ids)/4:
            experiment_outcomes = []
            for _ in range(100):
                random_sel = rd.sample(gene_list, size)
                x = 0
                for gene in random_sel:
                    if gene in myc_targets:
                        x += 1
                experiment_outcomes.append(x)
            p_values = []
            #calculate p values for each observed match
            for observed_match in observed_matches:
                p_value = sum(1 for experiment_match in experiment_outcomes if experiment_match >= observed_match)/ len(experiment_outcomes)
                p_values.append(p_value)
            above_baseline_count = sum(1 for p in p_values if p < 0.05)
            significant_motifs.append((binding_site, above_baseline_count))
    #save into csv file
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Motif", "Above Baseline Count"])
        csv_writer.writerows(significant_motifs)

def motif_plotter(gene_ids,size,bed_file,motif_list,fasta_file,stats=False,fixed_height=20,operator="OR"):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end),gene_name))
            gene_list.append(gene_name)
    search_range = 50 
    genes_with_motif = []

    for chromosome, start, end,gene_name in tss_data:
        sequence = genome[chromosome][start - search_range: end + search_range].seq
        if check_motifs_in_sequence(sequence, motif_list, operator):
            genes_with_motif.append(gene_name)
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end), gene_name))
            gene_list.append(gene_name)
    matches_percentages,observed_matches_percentages = [],[]
    for i in range(len(gene_ids) - size+1):
        genes_subset = gene_ids[i:i + size]
        matches_percentages.append(len(set(genes_subset) & set(genes_with_motif)))

    x = range(1, len(matches_percentages) + 1)
    observed_matches_percentages.extend(matches_percentages)  #Add observed matches for p-value calculation
    experiment_outcomes = []
    for _ in range(10000):  # Number of random samples (e.g., 10000)
        random_sel = rd.sample(gene_list, size)
        x = 0
        for gene in random_sel:
            if gene in genes_with_motif:
                x += 1
        experiment_outcomes.append(x)
    print("in whole genome:",len(genes_with_motif),round(len(genes_with_motif)/len(gene_list),2)*100,"percent")

    threshold = np.percentile(experiment_outcomes, 95)
    x_axis = range(1, len(matches_percentages) + 1)
    p_values = []
    for observed_match in observed_matches_percentages:
        p_value = sum(1 for experiment_match in experiment_outcomes if experiment_match >= observed_match) / len(experiment_outcomes)
        p_values.append(p_value)
    background_probability = {}
    background_probability_list = []
    for i in np.arange(-1,fixed_height):
        p_value = sum(1 for experiment_match in experiment_outcomes if experiment_match >= i) / len(experiment_outcomes)
        background_probability[i] = p_value
        background_probability_list.append(p_value)

    #X-axis
    plt.xlabel('Gene sets')
    plt.xticks(range(1, len(gene_ids), 10))
    plt.xlim(0, len(observed_matches_percentages))

    # Y-axis pvalue
    colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)]
    for i in range(len(background_probability_list)-1,-1,-1):
        if background_probability_list[i] > 0.05 and background_probability_list[i+1] < 0.05:
            signi_theshol = i
            break
    mid = ((signi_theshol+signi_theshol-1)/2/len(background_probability_list))
    positions = [0, mid, 1]

    cmap = LinearSegmentedColormap.from_list('custom_colormap', list(zip(positions, colors)))
    plt.ylim(-1, fixed_height)
    plt.fill_between(x_axis, -1, 0, facecolor='white')
    plt.yticks(range(-1, len(background_probability_list)-1),background_probability_list)
    gradient = np.linspace(1, 0, 256).reshape(-1, 1)
    gradient = np.hstack((gradient, gradient))

    # Remove the ylabel for P-value gradient
    plt.ylabel('')
    plt.imshow(gradient, extent=[0, len(observed_matches_percentages), 0, len(background_probability_list)], aspect='auto', cmap=cmap, alpha=0.75)

    # Y-axis nr of matches
    plt.twinx()
    plt.ylim(-1, fixed_height)
    plt.ylabel('Number of Matches')
    plt.plot(x_axis, observed_matches_percentages, linestyle='dotted', color='gray', linewidth=2, markersize=0)
    plt.plot(x_axis, observed_matches_percentages, marker='o', color="black", alpha=1, linestyle="", markersize=4)
    plt.yticks(range(-1, len(background_probability_list)-1),range(-1, len(background_probability_list)-1))
    plt.rcParams["figure.figsize"] = [16,9]

    #plot random distirbution statistics (median, standard deviation)
    if stats:
        mean = statistics.mean(experiment_outcomes)
        stdev = statistics.stdev(experiment_outcomes)
        plt.plot(x_axis, [mean] * len(x_axis), marker='o', color="grey", alpha=1, linestyle="", markersize=4)
        plt.errorbar(x_axis,[mean] * len(x_axis),yerr=stdev,color="grey",errorevery=1,fmt="none")
        plt.plot(x_axis, [mean+stdev] * len(x_axis), marker='o', color="#5A5A5A", alpha=1, linestyle="", markersize=4)
        plt.plot(x_axis, [mean-stdev] * len(x_axis), marker='o', color="#5A5A5A", alpha=1, linestyle="", markersize=4)
        plt.axhline(mean,linestyle="dashed",color="#5A5A5A")
        plt.axhline(mean+stdev,linestyle="dashdot",color="#5A5A5A")
        plt.axhline(mean-stdev,linestyle="dashdot",color="#5A5A5A")
    title =("Genes with %s motif"%(",".join(motif_list)))
    plt.title(title)
    plt.show()

def motif_plotter_old(gene_ids,size,bed_file,motif,fasta_file):
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Read the BED file
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            #print(line.strip().split("\t"))
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end),gene_name))
            gene_list.append(gene_name)
    search_range = 50 
    myc_targets = []

    for chromosome, start, end,gene_name in tss_data:
        sequence = genome[chromosome][start - search_range: end + search_range].seq
        if motif in sequence:
            myc_targets.append(gene_name)
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end), gene_name))
            gene_list.append(gene_name)
    genes_with_motif = myc_targets
    matches_percentages,observed_matches_percentages = [],[]
    for i in range(len(gene_ids) - size+1):
        genes_subset = gene_ids[i:i + size]
        matches_percentages.append(len(set(genes_subset) & set(genes_with_motif)))

    x = range(1, len(matches_percentages) + 1)
    plt.plot(x[:len(x) - 1], matches_percentages[:len(matches_percentages) - 1], marker='.')

    plt.xlabel('Gene sets')
    plt.ylabel('Number of Matches')
    plt.xticks(x)

    observed_matches_percentages.extend(matches_percentages)  # Add observed matches for p-value calculation

    experiment_outcomes = []
    for _ in range(10000):  # Number of random samples (e.g., 10000)
        random_sel = rd.sample(gene_list, size)
        x = 0
        for gene in random_sel:
            if gene in myc_targets:
                x += 1
        experiment_outcomes.append(x)

    plt.rcParams["figure.figsize"] = [16, 9]
    max_theshold = np.percentile(experiment_outcomes, 100)
    threshold = np.percentile(experiment_outcomes, 95)

    x_axis = range(1, len(matches_percentages) + 1)
    #plt.axhline(max_theshold, color="red", linestyle="dashed",mouseover=True)
    plt.fill_between(x_axis, threshold, 0, facecolor='red', alpha=0.5)

    # Calculate p-values for each observed match percentage
    p_values = []
    for observed_match in observed_matches_percentages:
        p_value = sum(1 for experiment_match in experiment_outcomes if experiment_match > observed_match) / len(experiment_outcomes)
        p_values.append(p_value)

    # Plot the observed matches and their p-values
    plt.plot(x_axis, observed_matches_percentages, marker='.', color="green")

    # Plot the p-values as a separate plot with a secondary y-axis
    plt.twinx()
    plt.plot(x_axis, p_values, color="blue",marker=".")
    plt.ylabel('P-values', color='blue')
    plt.tick_params(axis='y', colors='blue')
    plt.rcParams["figure.figsize"] = [16, 9]
    plt.xticks(range(0, len(gene_ids), 10))
    title = ("Genes in Choanocyte Trajectory with %s motif (pearson coeff 0.8)"% motif)
    plt.title(title)
    plt.show()

def target_genes(gene_ids,bed_file,motif,fasta_file,size=30,operator="OR"):
    significant_genes = []
    genes_subset_list = []
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # Read the BED file
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            #print(line.strip().split("\t"))
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end),gene_name))
            gene_list.append(gene_name)
    search_range = 50
    myc_targets = []

    for chromosome, start, end,gene_name in tss_data:
        sequence = genome[chromosome][start - search_range: end + search_range].seq
        if check_motifs_in_sequence(sequence, motif, operator):
            myc_targets.append(gene_name)
    tss_data = []
    gene_list = []
    with open(bed_file, "r") as file:
        for line in file:
            chromosome, start, end, gene_name = line.strip().split("\t")
            tss_data.append((chromosome, int(start), int(end), gene_name))
            gene_list.append(gene_name)
    genes_with_motif = myc_targets
    matches_percentages,observed_matches_percentages = [],[]
    for i in range(len(gene_ids) - size+1):
        genes_subset = gene_ids[i:i + size]
        matches_percentages.append(len(set(genes_subset) & set(genes_with_motif)))
        genes_subset_list.append(genes_subset)

    observed_matches_percentages.extend(matches_percentages)  # Add observed matches for p-value calculation

    experiment_outcomes = []
    for _ in range(10000):  # Number of random samples (e.g., 10000)
        random_sel = rd.sample(gene_list, size)
        x = 0
        for gene in random_sel:
            if gene in myc_targets:
                x += 1
        experiment_outcomes.append(x)

    plt.rcParams["figure.figsize"] = [16, 9]
    max_theshold = np.percentile(experiment_outcomes, 100)
    threshold = np.percentile(experiment_outcomes, 95)
    x_axis = range(1, len(matches_percentages) + 1)
    # Calculate p-values for each observed match percentage
    p_values = []
    for observed_match in observed_matches_percentages:
        p_value = sum(1 for experiment_match in experiment_outcomes if experiment_match > observed_match) / len(experiment_outcomes)
        p_values.append(p_value)

    for subset in genes_subset_list:
        if len(set(subset) & set(genes_with_motif)) >= threshold:
            for gene in subset:
                if gene in genes_with_motif:
                    significant_genes.append(gene)
    return(list(set(significant_genes)))