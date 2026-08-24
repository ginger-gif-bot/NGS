import os
from Bio import SeqIO
from pyliftover import LiftOver

# Analysing the reference file
ref = os.path.join("data","reference","e_coli_ref_genome","ncbi_dataset","data","GCF_000005845.2","genomic.gbff")

target_genes = ["rpoA","rpoB","rpoC"]
cds_list = []
for record in SeqIO.parse(ref,"genbank"):
    for feature in record.features:
        if feature.type == "CDS":
            gene = feature.qualifiers.get("gene",[""])[0]
            # if "rpo" in gene.lower():
            #     print(gene,feature.qualifiers.get("locus_tag",[""])[0]) #to find the correct gene name
            if gene in target_genes:
                cds_list.append({
                    "gene": gene,
                    "start":int(feature.location.start),
                    "end": int(feature.location.end),
                    "strand": int(feature.location.strand)
                })
print(cds_list)

# The liftover
chain_path = os.path.join("results","chain","SRR38870271.chain")
lo = LiftOver(chain_path)
gene = cds_list[0]
contig = "NC_000913.3"

start_lifted = lo.convert_coordinate(contig,gene["start"])
end_lifted = lo.convert_coordinate(contig,gene["end"])

print(start_lifted)
print(end_lifted)

# Slicing the consensus

consensus_path = os.path.join("results","consensus","SRR38870271.fasta")
consensus_record = SeqIO.read(consensus_path,"fasta")

start_pos = start_lifted[0][1]
end_pos = end_lifted[0][1]
cds_seq = consensus_record.seq[start_pos:end_pos]
if gene["strand"] == -1:
    cds_seq = cds_seq.reverse_complement()
print(f">{gene["gene"]}")
print(cds_seq[:60]) 