import argparse, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from read_genome import read_human_genome, read_input_genome
from download_genome import download_input_genomes
from results import add_to_list
from allignment_calc import process_contigs


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

    # Load human contigs
    human_contigs = list(read_human_genome(args.human))

    # Parallel processing
    with open(args.output, "w") as out_handle, ProcessPoolExecutor() as executor:
        for local_file in local_files:
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
