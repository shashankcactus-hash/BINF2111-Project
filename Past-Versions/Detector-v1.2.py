from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
import urllib.request, gzip, numpy as np
import argparse # New import for recomendations by Doc
#code changed to allow url vs url gzip and local files as well get the fna.gz URL from NCBI
# examples of running the script
#python Detector-v1.2.py Url.gz Url2.gz -o name_name2_fasta.md
#python Detector-1.2.py Local.fasta Url.gz -o local_url_fasta.md
percentage_list = []

# === percentage list and average ===
def add_to_list(percentage: float, out_handle, final=False):
    """Stores percentages, updates the average, and writes the results down."""
    
    # Explicitly declare intent to modify the global list 
    global percentage_list 
    
    if not final:
        percentage_list.append(round(percentage, 3)) 
        out_handle.write(f"This sample is likely {percentage:.2f}% the same!\n")  

        if percentage >= 95.00:
            out_handle.write("This sample is likely the same!\n \n")
        else:
            out_handle.write("This sample is likely not the same!\n \n")

    else:
        # Calculate average only if there is data to use!
        percent_average = np.mean(percentage_list) if percentage_list else 0.00
        out_handle.write(f"\nOn average, this sample is {percent_average:.2f}% the same!\n")


# === file/URL readers ===
def read_input_file(filename: str) -> SeqRecord:
    """Reads the local example FASTA file (must contain only one record)."""
    return SeqIO.read(filename, "fasta")

def read_fasta_source(source: str):
    """Reads a FASTA source, which can be a local file or a gzipped URL."""
    is_url = source.lower().startswith("http")
    
    if is_url:
        # Logic for reading from a URL
        try:
            with urllib.request.urlopen(source) as remote_file:
                # Check for .gz extension to decide on gzip processing
                if source.lower().endswith(".gz"):
                    with gzip.open(remote_file, "rt") as fasta:
                        # Use yield from to return a generator of SeqRecords
                        yield from SeqIO.parse(fasta, "fasta")
                else:
                    # For non-gzipped URLs
                    yield from SeqIO.parse(remote_file, "fasta")
        except Exception as e:
            print(f"Error reading FASTA from URL {source}: {e}")
            return
    else:
        # Logic for reading a local file (using SeqIO.parse for multi-record handling)
        try:
            yield from SeqIO.parse(source, "fasta")
        except FileNotFoundError:
            print(f"Error: Local file '{source}' not found.")
            return

# === alignment and percentage calculator ===
def alignment_calculator(output_seq: str, example_seq: str) -> float:
    """Performs local alignment and calculates percentage match."""
    
    segment_limit = 100000
    
    # The 'output' sequence (query) is limited to the first 100k bases (for time)
    query_segment = str(output_seq[:segment_limit])
    
    # The 'example' sequence (subject) is compared against this segment
    subject_segment = str(example_seq) # Use the whole subject sequence for comparison
    
    # Alignment 
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
                        
    score = aligner.score(query_segment, subject_segment)
    
    matches = score
    
    # The divisor should be the length of the segment being compared against 
    # (the limited 100k bases of the query sequence for saving TIME!)
    bases = len(query_segment) 
    
    if bases == 0:
        return 0.0
        
    percentage = (matches / bases) * 100
    return percentage

    
# === main === 
def main():
    # --- 1. Set up Argument Parsing ---
    parser = argparse.ArgumentParser(
        description="Compare the sequence similarity of a SUBJECT FASTA against a QUERY FASTA (file or URL)."
    )
    # The sequence being compared against (Reference/Query, e.g., Human genome)
    parser.add_argument(
        "query_source", 
        type=str, 
        help="The reference sequence source (e.g., human genome). Can be a local FASTA file or a remote URL (e.g., NCBI link)."
    )
    # The sequence being measured (Subject/Example, e.g., Chimp genome)
    parser.add_argument(
        "subject_source", 
        type=str, 
        help="The sequence to be measured/compared (e.g., animal genome). Can be a local FASTA file or a remote URL."
    )
    parser.add_argument(
        "-o", "--output", 
        type=str, 
        default="blastresult.md", 
        help="The output markdown file name. Default is blastresult.md."
    )
    
    args = parser.parse_args()
    
    querysource = args.query_source
    subjectsource = args.subject_source
    outputfile = args.output
    
    # --- 2. Process Subject Sequence (The sequence being measured) ---
    subject_records = list(read_fasta_source(subjectsource))
    if not subject_records:
        print(f"Could not read any records from subject source: {subjectsource}. Exiting.")
        return
        
    # We only use the first record of the subject for comparison
    subject = subject_records[0] 

    # --- 3. Process Query Sequence (The sequence being compared against) ---
    max_iterations = 6 # Limit iterations when querying large genomes

    with open(outputfile, "w") as out_handle:
        out_handle.write(f"# Comparison of SUBJECT: {subjectsource} (Record: {subject.id}) vs QUERY: {querysource}\n\n")

        # Iterate through the records of the QUERY source
        for i, query in enumerate(read_fasta_source(querysource)):
            
            if i >= max_iterations: 
                out_handle.write(f"\n--- Break: Limiting to {max_iterations} query records ---\n")
                break
                
            out_handle.write(f"{i+1}, Query ID: {query.id}, Length: {len(query.seq)}\n")

            # Calculates the percentage: (Query segment, Subject sequence)
            # Alignment is performed: Subject.seq is aligned against Query.seq[:100000]
            percentage = alignment_calculator(query.seq, subject.seq) 
            
            add_to_list(percentage, out_handle) 

        # Final average calculation and writing
        add_to_list(0, out_handle, final=True) 

# Execute the main function
if __name__ == "__main__":
    main()
