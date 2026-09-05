import os
from Bio import SeqIO
from glob import glob
from Bio.Data import CodonTable
from collections import Counter, defaultdict

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
gene_list = ["rpoB","katG","gyrA"]
valid_bases = set("ATCG")

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