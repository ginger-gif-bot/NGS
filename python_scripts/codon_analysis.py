# =============================================================================
# CODON ANALYSIS PIPELINE — codon_analysis.py
# Author: Kiran Gupta | M.Sc. Biotechnology Thesis
# Organism: Mycobacterium tuberculosis H37Rv (GCF_000195955.2)
# Dataset: CRyPTIC consortium — 50 drug-sensitive, 50 drug-resistant isolates
# Genes analysed: rpoB, katG, gyrA
#
# ⚠️  WARNING: This script reads 300 FASTA files and performs heavy computation.
#     Run ONCE to generate CSV outputs in results/.
#     Do NOT re-run just to regenerate plots — use analysis_notebook.ipynb instead.
#     Run plots ONE AT A TIME — running all plots together causes lag and
#     blank figure outputs due to matplotlib memory limitations on this machine.
#
# OUTPUT FILES (saved to results/CSVs):
#     enc_gc3_all_samples.csv  — ENC and GC3 per sample per gene
#     rscu_all_samples.csv     — Per-sample RSCU values (long format)
#     rscu_avg_rpoB.csv        — Average RSCU per group for rpoB
#     rscu_avg_katG.csv        — Average RSCU per group for katG
#     rscu_avg_gyrA.csv        — Average RSCU per group for gyrA
# =============================================================================

# --- SECTION 1: IMPORTS ---
import os
from Bio import SeqIO
from glob import glob
from Bio.Data import CodonTable
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

# --- SECTION 2: GLOBAL SETTINGS ---

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
gene_list = ["rpoB","katG","gyrA"]
valid_bases = set("ATCG")

aa_to_codons_std = defaultdict(list)
for codon, aa in codon_table.forward_table.items():
    aa_to_codons_std[aa].append(codon)
# print(aa_to_codons_std)

# --- SECTION 3: LOAD GROUP METADATA ---

with open(os.path.join("metadata","sensitive_ids.txt"),"r") as f:
    sensitive_ids = set(f.read().splitlines())

with open(os.path.join("metadata","resistant_ids.txt"),"r") as f:
    resistant_ids =  set(f.read().splitlines())

# --- SECTION 4: FUNCTION — codon_counts(gene, path_file) ---

def codon_counts(gene,path_file):
    codon_list = []
    skipped = 0
    record = SeqIO.read(path_file,"fasta")
    id = os.path.basename(path_file)
    id = id.replace(f"_{gene}","")
    for i in range(0,len(record),3):
        codon = str(record.seq[i:i+3])
        if len(codon) == 3 and all(base in valid_bases for base in codon):
            codon_list.append(codon)
        else:
            skipped +=1
    codon_dict = Counter(codon_list)
    return id, codon_dict, skipped

# --- SECTION 5: BUILD CODON COUNT DICTIONARY ---

all_codons = {}
all_skipped = {}
for gene in gene_list:
    path = glob(os.path.join("results","cds_seq",gene,"*.fasta"))
    for files in path:
        sample_id, codon_count,skips = codon_counts(gene,files)
        sample_id = sample_id.replace(".fasta","")
        if sample_id not in all_codons:
            all_codons[sample_id] = {}
        all_codons[sample_id][gene] = codon_count
        if sample_id not in all_skipped:
            all_skipped[sample_id] = {}
        all_skipped[sample_id][gene] = skips
    # print(codon_count)

# print(f"Total samples: {len(all_codons)}")
# print(f"Genes per sample: {list(all_codons[list(all_codons.keys())[0]].keys())}")
# print(all_skipped)

#  --- SECTION 6: FUNCTION — rscu(aa_dict_std, codon_count_dict) ---

#### === RSCU CALCULATION ===

def rscu(aa_dict_std,all_codon_dict):
    rscu_dict = {}
    for codon, count in all_codon_dict.items():
        aa = codon_table.forward_table.get(codon,None)
        # print(aa)
        if aa is None:
            continue
        synonymous_codons = aa_dict_std.get(aa,[])
        # print(syn_codons)
        total_amino_acid_count = sum(all_codon_dict.get(syn_codon,0) for syn_codon in synonymous_codons)
        # print(total_syn_count)
        rscu_val = round((count * len(synonymous_codons)) / (total_amino_acid_count),2) 
        rscu_dict[codon] = rscu_val
    return rscu_dict

# --- SECTION 7: BUILD RSCU DICTIONARY ---
all_codons_rscu = {}
for sample_id in all_codons:
    for cdn_count in all_codons[sample_id]:
        rscu_val = rscu(aa_to_codons_std,all_codons[sample_id][cdn_count])
        if sample_id not in all_codons_rscu:
            all_codons_rscu[sample_id] = {}
        all_codons_rscu[sample_id][cdn_count] = rscu_val

# print(all_codons_rscu)
first_sample = list(all_codons_rscu.keys())[0]
# print(f"Total samples: {len(all_codons_rscu)}")
# print(f"Genes in first sample: {list(all_codons_rscu[first_sample].keys())}")
# print(f"RSCU values for rpoB: {all_codons_rscu[first_sample]['rpoB']}")

# --- SECTION 8: ASSIGN GROUP LABELS ---
group_labels = {}

for sample_id in all_codons_rscu:
    if sample_id in sensitive_ids:
        group_labels[sample_id] = "sensitive"
    else:
        group_labels[sample_id] = "resistant"

# print(group_labels)
# print(sum(1 for v in group_labels.values() if v == "sensitive"))
# print(sum(1 for v in group_labels.values() if v == "resistant"))

# --- SECTION 9: AVERAGE RSCU PER GROUP ---

rscu_values_all = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))
rscu_avg = defaultdict(lambda:defaultdict(lambda:defaultdict(float)))

for sample_id in all_codons_rscu:
    for gene in all_codons_rscu[sample_id]:
        group = group_labels[sample_id]
        till_gene = all_codons_rscu[sample_id][gene]
        for codon in till_gene:
            rscu_values_all[group][gene][codon].append(till_gene[codon])
            
# print(rscu_values_all)

for grp in rscu_values_all:
    for gene in rscu_values_all[grp]:
        for codon in rscu_values_all[grp][gene]:
            rscu_values = rscu_values_all[grp][gene][codon]
            rscu_avg[grp][gene][codon] = round(sum(rscu_values)/len(rscu_values),2)

# print(rscu_avg["sensitive"]["rpoB"]["CTG"])
# print(rscu_avg["resistant"]["rpoB"]["CTG"])
# --- SECTION 10: BUILD FOLD GROUPS FOR ENC ---
#### === ENC CALCULATION ===
two_fold = []  
three_fold = []
four_fold = []
six_fold = []

for aa,codons in aa_to_codons_std.items():
    degenracy = len(codons)
    if degenracy == 1:
        continue
    elif degenracy == 2:
        two_fold.append(aa)
    elif degenracy == 3:
        three_fold.append(aa)
    elif degenracy == 4:
        four_fold.append(aa)
    elif degenracy == 6:
        six_fold.append(aa)   

# print(f"Two Fold: {len(two_fold)} -> {two_fold}")
# print(f"Three Fold: {len(three_fold)} -> {three_fold}")
# print(f"Four Fold: {len(four_fold)} -> {four_fold}")
# print(f"six Fold: {len(six_fold)} -> {six_fold}")

#  --- SECTION 11: FUNCTION — enc_calc(codon_counter, fold) ---

def enc_calc(codon_counter,fold):
    f_hat_values = []
    for aa in fold:
        aa_to_codons = aa_to_codons_std.get(aa,None)
        ni_val = []
        for codons in aa_to_codons:
            ni = codon_counter.get(codons,0)
            ni_val.append(ni)
        n = sum(ni_val)
        sum_sq = sum(x**2 for x in ni_val)
        if n <= 1:
            continue
        f = round((sum_sq -n ) / (n * (n -1)),2)
        f_hat_values.append(f)

    f_hat_values_avg = round(sum(x for x in f_hat_values)/len(f_hat_values),2)

    return f_hat_values_avg
# print(enc_calc(all_codons["ERR4810461"]["gyrA"],two_fold)) 
# print(aa_to_codons_std["F"])

enc_values = defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
for sample_id in group_labels:
    for gene in all_codons[sample_id]:
        f2 = enc_calc(all_codons[sample_id][gene],two_fold)
        f3 = enc_calc(all_codons[sample_id][gene],three_fold)
        f4 = enc_calc(all_codons[sample_id][gene],four_fold)
        f6 = enc_calc(all_codons[sample_id][gene],six_fold)

        if f2 and f3 and f4 and f6:
            enc = round(2 + (9/f2) + (1/f3) + (5/f4) + (3/f6),2)
        else:
            enc = None
        enc_values[group_labels[sample_id]][sample_id][gene] = enc
 

first_group = list(enc_values.keys())[0]
first_sample = list(enc_values[first_group].keys())[0]
# print(f"Groups: {list(enc_values.keys())}")
# print(f"First sample in {first_group}: {enc_values[first_group][first_sample]}") 

# --- SECTION 12: FUNCTION — gc3_calc(codon_counter) ---
#### === GC3 CALCULATION ===

def gc3_calc(codon_counter):
    gc = 0
    total = sum(codon_counter.values())
    for codon,count in codon_counter.items():
        if codon[2] == "G" or codon[2] == "C":
            gc += count

    return round(gc / total,4)

# print(gc3_calc(all_codons["ERR4810467"]["gyrA"]))

gc3_values = defaultdict(lambda:defaultdict(lambda:defaultdict(float)))
for sample_id in group_labels:
    for gene in all_codons[sample_id]:
        gc3_val = gc3_calc(all_codons[sample_id][gene])
        gc3_values[group_labels[sample_id]][sample_id][gene] = gc3_val

# print(gc3_values) 

# --- SECTION 13: BUILD ENC AND GC3 DICTIONARIES ---
#### === ENC - GC3 PLOT ===

# gyrA
enc_resistant_gyrA = []
enc_sensitive_gyrA = []
gc3_resistant_gyrA = []
gc3_sensitive_gyrA = []

# rpoB
enc_resistant_rpoB = []
enc_sensitive_rpoB = []
gc3_resistant_rpoB = []
gc3_sensitive_rpoB = []

# katG
enc_resistant_katG = []
enc_sensitive_katG = []
gc3_resistant_katG = []
gc3_sensitive_katG = []

# print(f"rpoB GC3 resistant sample: {gc3_resistant_rpoB[:5]}")
# print(f"rpoB GC3 sensitive sample: {gc3_sensitive_rpoB[:5]}")
# print(f"gyrA GC3 resistant sample: {gc3_resistant_gyrA[:5]}")

# print(f"rpoB enc resistant sample: {enc_resistant_rpoB[:5]}")
# print(f"rpoB enc sensitive sample: {enc_sensitive_rpoB[:5]}")
# print(f"gyrA enc resistant sample: {enc_resistant_gyrA[:5]}")

for sample_id in group_labels:
    grp = group_labels[sample_id]

    if "katG" in enc_values[grp][sample_id]:  
        if grp == "resistant":     
            enc_resistant_katG.append(enc_values[grp][sample_id]["katG"])
            gc3_resistant_katG.append(gc3_values[grp][sample_id]["katG"])
        else:
            enc_sensitive_katG.append(enc_values[grp][sample_id]["katG"])
            gc3_sensitive_katG.append(gc3_values[grp][sample_id]["katG"])

    if "gyrA" in enc_values[grp][sample_id]:
        if grp == "resistant":
            enc_resistant_gyrA.append(enc_values[grp][sample_id]["gyrA"])
            gc3_resistant_gyrA.append(gc3_values[grp][sample_id]["gyrA"])
        else:
            enc_sensitive_gyrA.append(enc_values[grp][sample_id]["gyrA"])
            gc3_sensitive_gyrA.append(gc3_values[grp][sample_id]["gyrA"])

    if "rpoB" in enc_values[grp][sample_id]:
        if grp == "resistant":
            enc_resistant_rpoB.append(enc_values[grp][sample_id]["rpoB"])
            gc3_resistant_rpoB.append(gc3_values[grp][sample_id]["rpoB"])
        else: 
            enc_sensitive_rpoB.append(enc_values[grp][sample_id]["rpoB"])
            gc3_sensitive_rpoB.append(gc3_values[grp][sample_id]["rpoB"])

# print(f"rpoB: {len(enc_resistant_rpoB)} resistant, {len(enc_sensitive_rpoB)} sensitive")
# print(f"katG: {len(enc_resistant_katG)} resistant, {len(enc_sensitive_katG)} sensitive")
# print(f"gyrA: {len(enc_resistant_gyrA)} resistant, {len(enc_sensitive_gyrA)} sensitive")

# --- SECTION 14: EXTRACT PLOT LISTS ---

# =============================================================================
# PLOTS — ⚠️  RUN ONE AT A TIME. Comment out others before running.
# Each plot block follows: create fig → draw → savefig → show → close
# =============================================================================

fig, axes  = plt.subplots(1,3,figsize=(15,5))

gc3_range = np.linspace(0.1,0.99,200)
enc_expected = 2 + gc3_range + (29 / (gc3_range**2 + (1 - gc3_range)**2))

plt_info = [("rpoB",enc_resistant_rpoB,enc_sensitive_rpoB,gc3_resistant_rpoB,gc3_sensitive_rpoB),
            ("katG",enc_resistant_katG,enc_sensitive_katG,gc3_resistant_katG,gc3_sensitive_katG),
            ("gyrA",enc_resistant_gyrA,enc_sensitive_gyrA,gc3_resistant_gyrA,gc3_sensitive_gyrA)]

# fig.suptitle("ENC vs GC3 - M. tuberculosis (rpoB,katG,gyrA)",fontsize=14,fontweight="bold")

# for i , (gene,enc_r,enc_s,gc3_r,gc3_s) in enumerate(plt_info):
#     ax = axes[i]
#     ax.plot(gc3_range,enc_expected,color="black",label="EXPECTED",linestyle="-",linewidth=1.5)
#     ax.scatter(gc3_s,enc_s,color="skyblue",label="Sensitive",alpha=0.4,s=80,marker="o",edgecolors="blue")
#     ax.scatter(gc3_r,enc_r,color="pink",label="Resistant",alpha=0.4,s=80,marker="^",edgecolors="maroon")
#     ax.set_title(gene,fontsize=15,fontweight="bold")
#     ax.set_xlabel("GC3",fontsize=11,fontweight="bold")
#     ax.set_ylabel("ENC",fontsize=11,fontweight="bold")
#     ax.set_ylim(20,61)
#     ax.set_xlim(0,1)
#     ax.grid(True,linestyle="--",alpha=0.4,color="#494a4a")
#     ax.legend(framealpha=0.7)
#     ax.set_facecolor("#e8f6f6")
# plt.tight_layout()
# plt.savefig(os.path.join("results","plots","enc_gc3_plot.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()


#### === ENC Boxplot ===

enc_info = [
    ("rpoB",enc_resistant_rpoB,enc_sensitive_rpoB),
    ("katG",enc_resistant_katG,enc_sensitive_katG),
    ("gyrA",enc_resistant_gyrA,enc_sensitive_gyrA)
]

# fig.suptitle("ENC Distribution-Sensitive vs Resistant M. tuberculosis",fontsize=14,fontweight="bold")
# for i, (gene,enc_r,enc_s) in enumerate(enc_info):
#     ax = axes[i]
#     ax.boxplot([enc_r,enc_s],
#                labels=["Resistant","Sensitive"],
#                patch_artist=True,
#                boxprops=dict(facecolor="#f181c0",alpha=0.7),
#                medianprops=dict(color="#6F2703",linewidth=2.5))
#     ax.set_title(gene,fontsize=15,fontweight="bold")
#     ax.set_xlabel("Groups",fontsize=11,fontweight="bold")
#     ax.set_ylabel("ENC",fontsize=11,fontweight="bold")
#     ax.grid(True,linestyle="--",alpha=0.4)
#     ax.set_facecolor("#e8f6f6")
#     ranges = {"rpoB": (33.8, 34.8), "katG": (39.4, 40.2), "gyrA": (37.8, 39.2)}
#     ax.set_ylim(ranges[gene])
# plt.tight_layout()
# plt.savefig(os.path.join("results","plots","enc_boxplot.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()

#### === GC3 Boxplot ===

gc3_info = [
    ("rpoB",gc3_resistant_rpoB,gc3_sensitive_rpoB),
    ("katG",gc3_resistant_katG,gc3_sensitive_katG),
    ("gyrA",gc3_resistant_gyrA,gc3_sensitive_gyrA)
]

# fig.suptitle("GC3 Distribution-Sensitive vs Resistant M. tuberculosis",fontsize=14,fontweight="bold")
# for i, (gene,gc3_r,gc3_s) in enumerate(gc3_info):
#     ax = axes[i]
#     ax.boxplot([gc3_r,gc3_s],
#                labels=["Resistant","Sensitive"],
#                patch_artist=True,
#                boxprops=dict(facecolor="#f181c0",alpha=0.7),
#                medianprops=dict(color="#6F2703",linewidth=2.5))
#     ax.set_title(gene,fontsize=15,fontweight="bold")
#     ax.set_xlabel("Groups",fontsize=11,fontweight="bold")
#     ax.set_ylabel("GC3",fontsize=11,fontweight="bold")
#     ax.grid(True,linestyle="--",alpha=0.4)
#     ax.set_facecolor("#e8f6f6")

# plt.tight_layout()
# plt.savefig(os.path.join("results","plots","gc3_boxplot.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()

#### === Heatmap for RSCU ===

data = {}
for codon in rscu_avg["resistant"]["rpoB"]:
    sens_val = rscu_avg["sensitive"]["rpoB"].get(codon,0)
    res_val = rscu_avg["resistant"]["rpoB"].get(codon,0)
    data[codon] = {"Sensitive":sens_val,"Resistant":res_val}

df_rpoB = pd.DataFrame(data).T 

for codon in rscu_avg["resistant"]["katG"]:
    sens_val = rscu_avg["sensitive"]["katG"].get(codon,0)
    res_val = rscu_avg["resistant"]["katG"].get(codon,0)
    data[codon] = {"Sensitive":sens_val,"Resistant":res_val}

df_katG = pd.DataFrame(data).T 

for codon in rscu_avg["resistant"]["gyrA"]:
    sens_val = rscu_avg["sensitive"]["gyrA"].get(codon,0)
    res_val = rscu_avg["resistant"]["gyrA"].get(codon,0)
    data[codon] = {"Sensitive":sens_val,"Resistant":res_val}

df_gyrA = pd.DataFrame(data).T 

fig,axes = plt.subplots(1,3,figsize=(10,16))
for i, (gene,df) in enumerate(zip(["rpoB","katG","gyrA"],[df_rpoB,df_katG,df_gyrA])):
    sns.heatmap(df,cmap="RdYlGn",annot=False,
            linewidths=0.5,ax=axes[i],vmin=0,vmax=4.5,linecolor="black")
    axes[i].set_title(gene,fontsize=13,fontweight="bold")
    axes[i].tick_params(axis="y",labelsize=7)
fig.suptitle("RSCU Heatmap - All Genes M. tuberculosis",fontsize=15,fontweight="bold")
plt.yticks(fontsize=6)
plt.tight_layout()
plt.savefig(os.path.join("results","plots","RSCU_all_genes_heatmap.png"),dpi=300,bbox_inches="tight")
plt.show()

plt.figure(figsize=(10,16))
# sns.heatmap(df_rpoB,cmap="RdYlGn",annot=False,
#             linewidths=1.5,vmin=0,vmax=4.5,linecolor="black")
# plt.title("RCSU Heatmap - rpoB",fontsize=13,fontweight="bold")
# plt.tight_layout()
# plt.yticks(fontsize=6)
# plt.savefig(os.path.join("results","plots","RSCU_ropB_heatmap.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()

# sns.heatmap(df_gyrA,cmap="RdYlGn",annot=False,
#             linewidths=1.5,vmin=0,vmax=4.5,linecolor="black")
# plt.title("RCSU Heatmap - gyrA",fontsize=13,fontweight="bold")
# plt.tight_layout()
# plt.yticks(fontsize=6)
# plt.savefig(os.path.join("results","plots","RSCU_gyrA_heatmap.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()

# sns.heatmap(df_katG,cmap="RdYlGn",annot=False,
#             linewidths=1.5,vmin=0,vmax=4.5,linecolor="black")
# plt.title("RCSU Heatmap - katG",fontsize=13,fontweight="bold")
# plt.tight_layout()
# plt.yticks(fontsize=6)
# plt.savefig(os.path.join("results","plots","RSCU_katG_heatmap.png"),dpi=300,bbox_inches="tight")
# plt.show()
# plt.close()

# --- SECTION 15: SAVE COMPUTED VALUES TO CSV ---

df_rpoB.to_csv(os.path.join("results","CSVs","rscu_avg_rpoB.csv"))
# print(df_rpoB.head())
df_katG.to_csv(os.path.join("results","CSVs","rscu_avg_katG.csv"))
# print(df_katG.head())
df_gyrA.to_csv(os.path.join("results","CSVs","rscu_avg_gyrA.csv"))
# print(df_gyrA.head())

rows = []
for sample_id in group_labels:
    grp = group_labels[sample_id]
    row = {"Sample_ID":sample_id,"Group":grp}
    for gene in gc3_values[grp][sample_id]:
        if gene in enc_values[grp][sample_id]:
            row[f"enc_{gene}"] = enc_values[grp][sample_id][gene]
            row[f"gc3_{gene}"] = gc3_values[grp][sample_id][gene] 
    rows.append(row)
df_results = pd.DataFrame(rows)
df_results.to_csv(os.path.join("results","CSVs","enc_gc3_all_samples.csv"),index=False)
# print(df_results.tail())

rscu_rows = []
for sample_id in all_codons_rscu:
    grp = group_labels[sample_id]
    for gene in all_codons_rscu[sample_id]:
        for codon,val in all_codons_rscu[sample_id][gene].items():
            rscu_rows.append({
                "Sample_Id":sample_id,
                "Group": grp,
                "Gene":gene,
                "Codon":codon,
                "RSCU":val
            })
df_rscu_all = pd.DataFrame(rscu_rows)
df_rscu_all.to_csv(os.path.join("results","CSVs","rscu_all_samples.csv"),index=False)
# print(df_rscu_all.head())
# print(df_rscu_all.shape)