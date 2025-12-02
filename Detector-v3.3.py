import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import urllib.request, gzip, os, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
import uuid
import urllib.parse

# === global percentage list ===
percentage_list = []

# === percentage list and average ===
def add_to_list(percentage, out_handle, final=False):
    if not final:
        percentage_list.append(round(percentage, 3))
        out_handle.write(f"This sample is {percentage:.2f}% human!\n")
        if percentage >= 95.0:
            out_handle.write("This sample is human!\n\n")
        else:
            out_handle.write("This sample is not human!\n\n")
    else:
        if percentage_list:
            percent_average = np.mean(percentage_list)
            out_handle.write(f"\nOn average, this sample is {percent_average:.2f}% human!\n")
        else:
            out_handle.write("\nNo valid alignments were performed.\n")

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
            for i, record in enumerate(subset[:max_sequences]):
                yield record

# === alignment calculator ===
def alignment_calculator(human_seq, input_seq, max_bases=100000):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    human_slice = str(human_seq[:max_bases])
    input_slice = str(input_seq[:max_bases])
    score = aligner.score(human_slice, input_slice)
    percentage = (score / len(input_slice)) * 100
    return percentage

# === process one input contig against human contigs ===
def process_contigs(input_seq, human_contigs):
    local_results = []
    for human_seq in human_contigs:
        pct = alignment_calculator(human_seq, input_seq)
        local_results.append((input_seq.id, human_seq.id, pct))
    return local_results

# === download input genome(s) with validation ===
def download_input_genomes(urls):
    """
    Downloads a list of genome URLs ending with '.fna.gz'.
    Saves them locally using the last part of the URL (e.g., GCA_*.fna.gz)
    Skips download if file exists.
    Returns list of local filenames.
    """
    local_files = []
    for url in urls:
        if not url.endswith(".fna.gz"):
            raise ValueError(f"URL must end with .fna.gz: {url}")

        # Use the filename from the URL itself
        local_file = os.path.basename(urllib.parse.urlparse(url).path)

        if os.path.exists(local_file):
            print(f"{local_file} already exists, skipping download.")
        else:
            print(f"Downloading {url} -> {local_file}")
            urllib.request.urlretrieve(url, local_file)
            print("Download complete!")

        local_files.append(local_file)

    return local_files

# === main ===
def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Compare genome(s) to human genome.")

    parser.add_argument(
        "inputseqs",
        nargs="+",
        help="One or more URLs pointing to gzipped FASTA genomes (.fna.gz)."
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        default="blastresult.md",
        help="Output markdown file name."
    )

    parser.add_argument(
        "--human",
        type=str,
        default="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        help="Human genome reference (default GRCh38.p14)."
    )

    args = parser.parse_args()

    # Download & validate input genomes
    local_files = download_input_genomes(args.inputseqs)

    # Load human contigs once
    human_contigs = list(read_human_genome(args.human))

    # Parallel processing
    with open(args.output, "w") as out_handle, ProcessPoolExecutor() as executor:
        for local_file in local_files:
            percentage_list.clear()  # reset percentages per genome
            out_handle.write(f"\n=== Processing {local_file} ===\n")
            futures = {executor.submit(process_contigs, input_seq, human_contigs): input_seq.id
                       for input_seq in read_input_genome(local_file, max_contigs=5)}

            for future in as_completed(futures):
                input_id_placeholder = futures[future]
                try:
                    process_results = future.result()
                    if process_results:
                        current_input_id = process_results[0][0]
                        out_handle.write(f"\nInput contig {current_input_id} results\n")
                        for input_id, human_id, pct in process_results:
                            out_handle.write(f"Against Human contig {human_id}: \n")
                            add_to_list(pct, out_handle)
                    else:
                        out_handle.write(f"\nNo results for input contig {input_id_placeholder}.\n")
                except Exception as e:
                    out_handle.write(f"Error processing {input_id_placeholder}: {str(e)}\n")

            add_to_list(0, out_handle, final=True)

    end_time = time.time()
    print(f"The total runtime is {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()
