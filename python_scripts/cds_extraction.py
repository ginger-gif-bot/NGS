import os
from Bio import SeqIO
from pyliftover import LiftOver
from glob import glob
from datetime import datetime
from time import time

ref = os.path.join("data","reference","tb_ref_genome","ncbi_dataset","data","GCF_000195955.2","genomic.gbff")

def gene_name_extraction(ref,gene_name):
    gene_names_list = []
    for records in SeqIO.parse(ref,"genbank"):
        for feature in records.features:
            if feature.type == "CDS":
                gene = feature.qualifiers.get("gene",[""])[0]
                if gene_name in gene.lower():
                    gene_names_list.append({
                        "gene_name": gene,
                        "locus_tag": feature.qualifiers.get("locus_tag",[""])[0]
                    })
    return gene_names_list

# gene_names for TB = rpoB,katG,gyrA
gene_list = gene_name_extraction(ref,"gyra")
# print(gene_list)
def cds_list_extraction(gene_list):
    cds_list = []
    for gene in gene_list:
        cds_list.append(gene["gene_name"])
    return cds_list

cds_list = cds_list_extraction(gene_list)
target_list = ["rpoB","katG","gyrA"]

def cds_info(target_gene_list,ref):
    cds_info = []
    for records in SeqIO.parse(ref,"genbank"):
        for feature in records.features:
            if feature.type == "CDS":
                gene = feature.qualifiers.get("gene",[""])[0]
                if gene in target_gene_list:
                    cds_info.append({
                        "gene":gene,
                        "contig": records.id,
                        "start":int(feature.location.start),
                        "end": int(feature.location.end),
                        "strand": int(feature.location.strand)
                    })
    
    return cds_info

cds_infomation = cds_info(target_list,ref) 
# print(cds_infomation)
log_file = open(os.path.join("logs","cds_extraction.log"),"a")

def log(msg):
    print(msg)
    log_file.write(msg + "\n")
    log_file.flush()

log(f"=== Run started {datetime.now()} ===")

def cds_extraction(cds_info):
    os.makedirs(os.path.join("results","cds_seq"),exist_ok=True)
    

    for consensus_path in glob(os.path.join("results","consensus","*.fasta")):
        sample_id = os.path.basename(consensus_path).replace(".fasta","")
        chain_path = os.path.join("results","chain",f"{sample_id}.chain")

        lo = LiftOver(chain_path)
        consensus_record = SeqIO.read(consensus_path,"fasta")

        for sample in cds_info:
            gene = sample["gene"]
            gene_dir = os.path.join("results","cds_seq",f"{gene}")
            os.makedirs(gene_dir,exist_ok=True)
            out_path = os.path.join(gene_dir,f"{sample_id}_{gene}.fasta")
            
            if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
                log(f"Skipping {sample_id}_{gene}, already exists")
                continue

            start_time = time()
            contig = sample["contig"]
            start_lifted = lo.convert_coordinate(contig,sample["start"])
            end_lifted = lo.convert_coordinate(contig,sample["end"])
            start_pos = start_lifted[0][1]
            end_pos = end_lifted[0][1]
            cds_seq = consensus_record.seq[start_pos:end_pos]
            
            if sample["strand"] == -1:
                cds_seq = cds_seq.reverse_complement()
 
            with open(out_path,"w") as f:
                f.write(f">{sample_id}_{gene}\n{str(cds_seq)}\n")

            elapsed = time() - start_time
            seq_len = end_pos - start_pos
            log(f"{sample_id}_{gene} done in {elapsed:.2f}s ({seq_len%3} frame check)")
            

print(cds_extraction(cds_infomation))

def file_check(target_gene):
    files_in_folder = []
    for gene in target_gene:
        files = glob(os.path.join("results","cds_seq",gene,"*.fasta"))
        files_in_folder.append({
            gene : len(files)
        })
    return files_in_folder

print(file_check(target_list))

log_file.close() 