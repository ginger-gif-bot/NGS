import os
from Bio import SeqIO
from glob import glob

gene_list = ["rpoB","katG","gyrA"]

count = 0
excluded = {}
for gene in gene_list:
    path = glob(os.path.join("results","cds_seq",gene,"*.fasta"))
    for files in path:
        # print(files)
        record = SeqIO.read(files,"fasta")
        if len(record.seq)%3 != 0:
            print(f"Frame Issue: {files}")
            if gene not in excluded:
                excluded[gene] = []
            excluded[gene].append(files)
        else:
            count +=1
print(f"{count} files OK")

print(excluded)

# To move the files with sequence length not divisible by 3 into another folder

for gene,file_list in excluded.items():
    excluded_dir = os.path.join("results","cds_seq",gene,f"excluded_{gene}")
    os.makedirs(excluded_dir,exist_ok=True)
    for files in file_list:
        filename = os.path.basename(files)
        destination = os.path.join(excluded_dir,filename)
        os.rename(files,destination)
        print(f"Moved {filename} -> {excluded_dir}")


# To check for duplicate files
id_list = {}
for gene in gene_list:
    path = glob(os.path.join("results","cds_seq",gene,"*.fasta"))
    for files in path:
        id = os.path.basename(files)
        id = id.replace(f"_{gene}","")
        if gene not in id_list:
            id_list[gene] = []
        id_list[gene].append(id)

        if len(id_list[gene]) == len(set(id_list[gene])):
            print("No Duplicates found")
        else:
            print(f"{id} duplicte found")

# print(id_list["gyrA"])