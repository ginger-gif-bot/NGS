import os
from Bio import SeqIO
from glob import glob
from Bio.Data import CodonTable
from collections import Counter, defaultdict

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
gene_list = ["rpoB","katG","gyrA"]
valid_bases = set("ATCG")

with open(os.path.join("metadata","sensitive_ids.txt"),"r") as f:
    sensitive_ids = set(f.read().splitlines())

with open(os.path.join("metadata","resistant_ids.txt"),"r") as f:
    resistant_ids =  set(f.read().splitlines())

#### === RSCU CALCULATION ===

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

print(f"Total samples: {len(all_codons)}")
print(f"Genes per sample: {list(all_codons[list(all_codons.keys())[0]].keys())}")
# print(all_skipped)

aa_to_codons_std = defaultdict(list)
for codon, aa in codon_table.forward_table.items():
    aa_to_codons_std[aa].append(codon)

# print(aa_to_codons_std)

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

all_codons_rscu = {}
for sample_id in all_codons:
    for cdn_count in all_codons[sample_id]:
        rscu_val = rscu(aa_to_codons_std,all_codons[sample_id][cdn_count])
        if sample_id not in all_codons_rscu:
            all_codons_rscu[sample_id] = {}
        all_codons_rscu[sample_id][cdn_count] = rscu_val

# print(all_codons_rscu)
first_sample = list(all_codons_rscu.keys())[0]
print(f"Total samples: {len(all_codons_rscu)}")
print(f"Genes in first sample: {list(all_codons_rscu[first_sample].keys())}")
print(f"RSCU values for rpoB: {all_codons_rscu[first_sample]['rpoB']}")

group_labels = {}

for sample_id in all_codons_rscu:
    if sample_id in sensitive_ids:
        group_labels[sample_id] = "sensitive"
    else:
        group_labels[sample_id] = "resistant"

# print(group_labels)
print(sum(1 for v in group_labels.values() if v == "sensitive"))
print(sum(1 for v in group_labels.values() if v == "resistant"))

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

print(rscu_avg["sensitive"]["rpoB"]["CTG"])
print(rscu_avg["resistant"]["rpoB"]["CTG"])

#### === ENC CALCULATION
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

print(f"Two Fold: {len(two_fold)} -> {two_fold}")
print(f"Three Fold: {len(three_fold)} -> {three_fold}")
print(f"Four Fold: {len(four_fold)} -> {four_fold}")
print(f"six Fold: {len(six_fold)} -> {six_fold}")

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
print(f"Groups: {list(enc_values.keys())}")
print(f"First sample in {first_group}: {enc_values[first_group][first_sample]}") 

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