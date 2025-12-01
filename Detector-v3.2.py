import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import urllib.request, gzip, os, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess # new import for v3.0 for running Minimap2
import re # for parsing new strings
import uuid # for unique temp file names


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
            if i >= max_contigs:  # stop after first 2 contigs
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


# === process one input contig against an x amount of human contigs ===
def process_contigs(input_seq, human_contigs):
    # local list to avoid cross-process issues
    local_results = []
    for human_seq in human_contigs:
        # 🟢 FIX: Pass the full SeqRecord objects, not just the .seq attribute
        pct = alignment_calculator(human_seq, input_seq)
        local_results.append((input_seq.id, human_seq.id, pct))
    return local_results #returning results only calculated by this process

# === download input genome if it is missing ===
def ensure_input_genome_exists(url, local_file):
    if not os.path.exists(local_file):
        print(f"Downloading the input genome to {local_file} (this may take a while...)")
    else:
        print("Input genome already exists locally, and is saved!")
    urllib.request.urlretrieve(url, local_file)
    print("Download complete!")


# === main ===
def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Program to compare the sequence similarity of a given fasta with the human genome.")
    
    parser.add_argument(
        "inputseq",
        type=str,
        help="The reference sequence source (e.g., human genome). Can be a local FASTA file or a remote URL (e.g., NCBI link)."
    )
    
    parser.add_argument(
        "-i", "--inputfile",
        type=str,
        default="input.fna.gz",
        help="Filename to save the input sequence as. Default is input.fna.gz"
    )
    
    parser.add_argument(
        "-o", "--output", 
        type=str, 
        default="blastresult.md", 
        help="The output markdown file name. Default is blastresult.md."
    )

    parser.add_argument(
        "--human",
        type=str,
        default="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        help="A reference to the URL containing the human genome. default is GRCh38.p14"
    )
    args = parser.parse_args()

    inputseq = args.inputseq
    inputfile = args.inputfile
    humanseq = args.human
    outputfile = args.output
    #Ensure input genome exists locally
    ensure_input_genome_exists(inputseq, inputfile)

    #Read human contigs once
    human_contigs = list(read_human_genome(humanseq))

    #Parallel processing some changes to fit v3
    with open(outputfile, "w") as out_handle, ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_contigs, input_seq, human_contigs): input_seq.id
                   for input_seq in read_input_genome(inputfile, max_contigs=5)}

        for future in as_completed(futures):
            # 🟢 FIX: Define the tracking ID outside the try block for error handling
            input_id_placeholder = futures[future]
            
            try:
                process_results = future.result()
                
                # Check if process_results is valid before accessing elements
                if process_results:
                    current_input_id = process_results[0][0]
                    out_handle.write(f"\nInput contig {current_input_id} results ")
                    for input_id, human_id, pct in process_results:
                        out_handle.write(f"against Human contig {human_id}: \n")
                        #modifies the local percentage list
                        add_to_list(pct, out_handle)
                else:
                    out_handle.write(f"\nNo results for input contig {input_id_placeholder}.\n")
                    
            except Exception as e:
                # 🟢 FIX: Use the defined input_id_placeholder in the error message
                out_handle.write(f"Error processing {input_id_placeholder}: {str(e)}\n")

        add_to_list(0, out_handle, final=True)

        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f"The total runtime is {elapsed_time:.2f} seconds.")
    
if __name__ == "__main__":
    main()
