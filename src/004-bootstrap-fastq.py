"""
2020-08-13
GC

Boostrap resample a fastq file. Seqtk cannot be used since it performs sampling without
replacement and we need sampling with replacement. The script needs uncompressed paired-end
fastq files as input. It assumes that the identifiers used in R1 match those in R2. 

This program was written somewhat with memory usage in mind. Instead of adding every new reocrd
for a given replicate to an in-memory list and then writing the list out all at once, which is
potentially faster, the program indexes the fastq files and then iterates through the ids 
pulling out each record one at a time. 

TODO: 

- Modify for single end reads 
"""
from random import choices, seed
from pathlib import Path
import argparse

from Bio import SeqIO


parser = argparse.ArgumentParser(description="Perform boostrap resampling of fastq file")
parser.add_argument('fq1', help="Path to read 1 in uncompressed fastq format")
parser.add_argument('fq2', help="Path to read 2 in uncompressed fastq format")
parser.add_argument("--seed", help="Optional random seed. Defaults to system time", default=None)
parser.add_argument('--outDir', help="Optional directory to save the boostrapped fastq files", default=".")
parser.add_argument('--replicates', help="Optional number of bootstrap replicates to produce", default=1)
parser.add_argument('--k', help="Optional number of seqs to sample with replacement. Defaults to all seqs in fastq file")
args = parser.parse_args()

fq1_path = args.fq1
fq2_path = args.fq2
replicates = int(args.replicates)
out_dir = args.outDir
rand_seed = args.seed
k = args.k

if out_dir:
    out_dir = Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir()
else:
    out_dir = Path(".")

R1_name = Path(fq1_path).stem
R2_name = Path(fq2_path).stem

# main boostrap script --------------------------------------------------------
print(f"Indexing {fq1_path}...")
fq1 = SeqIO.index(fq1_path, "fastq")

print(f"Indexing {fq2_path}...")
fq2 = SeqIO.index(fq2_path, "fastq")

print(f"Generating bootstrap sample of record IDs...")
seed(rand_seed)
seq_ids = list(fq1.keys())

K = len(seq_ids) if k is None else int(k)
boot_ids = [choices(seq_ids, k = K) for r in range(replicates)]


def write_records(boot_list, fq_name, fq_dict):
    """Iterate over the bootstrap ids, writing a new fastq file for each
    replicate in the boot list.
    """
    for i, boots in enumerate(boot_list, start=1):
        fname = f"{fq_name}.boot.{str(i)}.fastq"
        outfile = Path(out_dir, fname)
        print(f"Writing boostrapped fastq file to {outfile}...")
        with open(outfile, "w") as out:
            for boot in boots:
                record = fq_dict[boot].format("fastq")
                out.write(record)

write_records(boot_list=boot_ids, fq_name=R1_name, fq_dict=fq1)
write_records(boot_list=boot_ids, fq_name=R2_name, fq_dict=fq2)
print("Done.")
