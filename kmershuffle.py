from collections import defaultdict
import pysam 
import pybedtools
import numpy as np
import pandas as pd
import argparse
import tempfile
import os

parser = argparse.ArgumentParser(description='Z-scoring of shuffled k-mers')

parser.add_argument('--anno', action="store", dest='anno', required=True)
parser.add_argument('--binds', action="store", dest='binds', required=True)
parser.add_argument('--seq', action="store", dest='seq', required=True)
parser.add_argument('--k', action="store", dest='k', default=3)
parser.add_argument('--shufflings', action="store", dest='shufflings', default=100)

args = parser.parse_args()
gff = args.anno
if not os.path.exists(gff):
    raise FileNotFoundError("GTF-file does not exist")
bed = args.binds
if not os.path.exists(bed):
    raise FileNotFoundError("BED-file does not exist")
fasta = args.seq
if not os.path.exists(fasta):
    raise FileNotFoundError("FASTA-file does not exist")
k = int(args.k)
shufflings = int(args.shufflings)


def count_kmers(seqs, k):
     kmers = {}
     sum_seqs = 0
     for seq in seqs:
        # rna = seq.replace("T", "U").upper().strip()
        for i in range(len(seq)-k+1):
            if seq[i:k] not in kmers.keys():
                kmers[seq[i:i+k]] = 0

            kmers[seq[i:i+k]] += 1
     return kmers

df = pd.read_table(gff, names=["chr",2,'feature',"start",'stop',6,"strand",8,'attr'])
new = df['attr'].str.split("\"", n=2, expand=True)
df["gene"] = new[1]
df.drop(columns=['attr'], inplace=True)
genes = df.groupby(["gene", 'strand'],as_index=False).agg({
    'chr': 'min',
    'start': 'min',
    'stop': 'max'
})
genes = genes[genes["gene"].str.fullmatch(r'(AT\dG\d{5})')==True].iloc[:, [2,3,4,1,0]]

site_length = 0

bind_sites = pybedtools.BedTool(bed)
bind_sites_small = []
dd = defaultdict(list)
kmers = {}
seqs = bind_sites.sequence(fi=fasta)
with tempfile.TemporaryDirectory() as tmpdir:
    seqs.save_seqs(f"{tmpdir}/binding_features.fa")
    with pysam.FastxFile(f"{tmpdir}/binding_features.fa") as fh:
        site_length = len(fh.__next__().sequence)
        kmers = count_kmers([f.sequence for f in fh],k)

    for feature in bind_sites:
        median = int(np.median(np.arange(feature.start, feature.stop)))
        bind_sites_small.append((feature.chrom, median, median))

    bind_features = pybedtools.BedTool.from_dataframe(genes[genes["gene"].str.fullmatch(r'(AT\dG\d{5})')==True]).intersect(bind_sites_small,u=True)
    count_features = bind_features.intersect(bind_sites_small, c=True)

    for shuffle in range(shufflings):
        print(f"Shuffle: {shuffle+1}")
        seqs = pybedtools.BedTool(list(np.concatenate([[(feature.chrom, site_start, site_start+site_length) for site_start in np.random.randint(feature.start, feature.stop+1, size=feature.count)] for feature in count_features]))).sequence(fi="TAIR10_chr_all.fa")
        seqs.save_seqs(f"{tmpdir}/binding_features{shuffle}.fa", )
        with pysam.FastxFile(f"{tmpdir}/binding_features{shuffle}.fa") as fh:
                shuffled_kmers = count_kmers([f.sequence for f in fh],k)
                for key, value in shuffled_kmers.items():
                    dd[key].append(value)

kmer_stats = {}
for key in kmers.keys():
    kmer_stats[key.replace("T","U")] = [kmers[key],np.mean(dd[key]), np.std(dd[key]), (kmers[key]-np.mean(dd[key]))/np.std(dd[key])]

stats = pd.DataFrame.from_dict(kmer_stats, orient="index", columns=["count","mean shufflings", "std shufflings", "z-score"])
stats.to_csv("kmer_stats.csv")
