from Bio import SeqIO
import urllib.request, gzip

# === read input genome locally ===
def read_input_genome(filename, max_contigs=2):
    with gzip.open(filename, "rt") as handle:
        for i, seq_record in enumerate(SeqIO.parse(handle, "fasta")):
            if i >= max_contigs:
                break
            yield seq_record

# === human genome streaming from NCBI ===
def read_human_genome(human_url, limit_fraction=0.2, max_sequences=3):
    with urllib.request.urlopen(human_url) as response:
        with gzip.open(response, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            subset_size = max(1, int(len(records) * limit_fraction))
            subset = records[:subset_size]
            for record in subset[:max_sequences]:
                yield record
