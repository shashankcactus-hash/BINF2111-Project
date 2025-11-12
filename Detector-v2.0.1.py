from Bio import SeqIO
from Bio.Align import PairwiseAligner
import urllib.request, gzip, os, time, numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import re


percentage_list = []
results = []

# === percentage list and average ===
def add_to_list(percentage, out_handle, final=False):
    if not final:
        percentage_list.append(round(percentage, 3))
        out_handle.write(f"This animal is {percentage:.2f}% human!\n")
        if percentage >= 95.0:
            out_handle.write("This animal is human!\n\n")
        else:
            out_handle.write("This animal is not human!\n\n")
    else:
        if percentage_list:
            percent_average = np.mean(percentage_list)
            out_handle.write(f"\nOn average, this animal is {percent_average:.2f}% human!\n")
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
#*this part should be a little quicker now... not by much i'd imagine
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
    for human_seq in human_contigs:
        pct = alignment_calculator(human_seq.seq, input_seq.seq)
        results.append((input_seq.id, human_seq.id, pct))
    return results

# === download input genome if it is missing ===
def ensure_input_genome_exists(url, local_file):
    if not os.path.exists(local_file):
        print(f"Downloading the input genome to {local_file} (this may take a while...)")
        urllib.request.urlretrieve(url, local_file)
        print("Download complete!")
    else:
        print("Input genome already exists locally, and is saved!")

# === displays runtime === 
#*please note this doesn't work as intended 
#*will display 0 seconds of runtime, but true runtime is around 2:35 


# === main ===
def main():
    starttime = time.time()

    inputseq = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/350/195/GCF_018350195.1_P.tigris_Pti1_mat1.1/GCF_018350195.1_P.tigris_Pti1_mat1.1_genomic.fna.gz"

    for i in range(1, 10000):
        inputfile = f"input{i}.fna.gz"
        if not os.path.exists(inputfile):
            print(f"Saving new file as {inputfile}!")
            break #i know not all inputs are going to be a .gz so if we can improve on how to fix it that would be great

    humanseq = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    outputfile = "blastresult.md"

    #Ensure input genome exists locally
    ensure_input_genome_exists(inputseq, inputfile)

    #Read human contigs once
    human_contigs = list(read_human_genome(humanseq))

    #Parallel processing
    with open(outputfile, "w") as out_handle, ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_contigs, input_seq, human_contigs): input_seq.id
                   for input_seq in read_input_genome(inputfile, max_contigs=5)}

        for future in as_completed(futures):
            input_id = futures[future]
            try:
                results = future.result()
                out_handle.write(f"\nInput contig {input_id} results ")
                for input_id, human_id, pct in results:
                    out_handle.write(f"against Human contig {human_id}: \n")
                    add_to_list(pct, out_handle)
            except Exception as e:
                out_handle.write(f"Error processing {input_id}: {str(e)}\n")

        add_to_list(0, out_handle, final=True)

    endtime = time.time()

    elapsedtime = endtime - starttime

    print(f"The total runtime is {elapsedtime:.2f} seconds.")
    
if __name__ == "__main__":
    main()
