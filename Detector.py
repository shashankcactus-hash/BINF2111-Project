from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import urllib.request, time, gzip
from io import StringIO



def human_percentage_blast():
    examplefasta = "chimp.fasta"  #File that is being used to compare to human. Currently its being used to compare and search for itself.

    humanfasta = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" #Human genome 3.1GB file
    
    outputfile = "blastresult.md" #All the results will get thrown into here, is markdown(md); light weight.

    percentage_list = []

    example = SeqIO.read(examplefasta, "fasta")

    with open(outputfile, "w") as out_handle: #As you are writing in the outputfile:
        with urllib.request.urlopen(humanfasta) as huma:
            with gzip.open(huma, "rt") as human:
                for i, hum in enumerate(SeqIO.parse(human, "fasta")):
                    out_handle.write(f"{i+1}, ID: {hum.id}, SeqRecord: {len(hum.seq)}\n")


                    #Alignment shit
                    aligner = PairwiseAligner()
                    aligner.mode = "local"
                    aligner.match_score = 1
                    aligner.mismatch_score = 0
                    aligner.open_gap_score = 0
                    aligner.extend_gap_score = 0


                    #Score calculator
                    matches = 0

                    bases = 0

                    score = aligner.score(str(hum.seq[:10000]), str(example.seq))

                    matches += score

                    bases += len(example.seq)

                    percentage = (matches / bases) * 100


                    #Stores all the values in a list

                    percentage_list.append(round(percentage, 3))


                    #Writes the shit into the md file
                    out_handle.write(f"This sample is {percentage:.2f}% human!\n")

                    if percentage != 100.00:
                        out_handle.write(f"This sample is not human!\n \n")
                    else:
                        out_handle.write(f"This sample is human!\n \n")


                    #So we don't brick our pc :D
                    if i >= 4:
                        break

human_percentage_blast() #Prints the class


    





