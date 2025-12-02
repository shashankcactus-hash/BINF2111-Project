from Bio import SeqIO
from Bio.Align import PairwiseAligner
import urllib.request, gzip, numpy as np


percentage_list = []


#Stores percentages, updates the average, and writes the results down in the md file

# === percentage list and average ===
def add_to_list(percentage, out_handle, final = False):

    if not final:
        percentage_list.append(round(percentage, 3)) #Rounds to 3 decimals places
        out_handle.write(f"This animal is {percentage:.2f}% human!\n") 

        if percentage >= 95.00:
            out_handle.write(f"This animal is human!\n \n")

        else:
            out_handle.write(f"This animal is not human!\n \n")

    else:
        percent_average = np.mean(percentage_list)
        out_handle.write(f"\nOn average, this animal is {percent_average:.2f}% human!\n") #Why is this guy not writing :rage:


# === user file ===
def read_input_file(filename):
    return SeqIO.read(filename, "fasta")


# === human genome file ===
def read_human_fasta(humanfasta):
        with urllib.request.urlopen(humanfasta ) as human:
            with gzip.open(human, "rt") as fasta:
                 for humanF in SeqIO.parse(fasta, "fasta"):
                     yield humanF


# === alignment and percentage calculator ===
def alignment_calculator(output_seq, example_seq):
    #Alignment 
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
                    
                    
    #Score calculator 
    matches = 0
    bases = 0
    score = aligner.score(str(output_seq[:100000]), str(example_seq))
    matches += score
    bases += len(example_seq)
    percentage = (matches / bases) * 100
    return percentage


    
# === main === 
def main():
    examplefasta = "chimp.fasta" #File that is being used to compare to human. Currently its being used to compare and search for itself.

    humanfasta = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" #Human genome 3.1GB file
    
    outputfile = "blastresult.md" #All the results will get thrown into here, is markdown(md); light weight.

    example = read_input_file(examplefasta) #Reads the example fasta file

    with open(outputfile, "w") as out_handle:
        for i, output in enumerate(read_human_fasta(humanfasta)):
            out_handle.write(f"{i+1}, ID: {output.id}, SeqRecord: {len(output.seq)}\n")

            percentage = alignment_calculator(output.seq, example.seq) #Calculates the percentage of human DNA in the example file
            add_to_list(percentage, out_handle) #Adds the percentage to the list and writes it to the md file


            #So we don't brick our pc :D
            if i >= 5: #For some reason always does 2 more than this number, I don't know why. 
                break
            
        add_to_list(0, out_handle, final = True) #Final average calculation             

main() #Prints the class

if __name__ == "__main__":
    main()
