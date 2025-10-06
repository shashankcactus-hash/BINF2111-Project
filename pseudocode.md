Import biopython methods  

initialize an empty list, percentage_list.  

define the function "add_to_list"  

if the iteration is not the final iteration, then proceed
add the percentage rounded to 3 decimals to percentage_list.
write this percentage out to a provided file.
checks if the percentage value is greater than or equal to 95%  

if true, then write to the provided file stating that the sample is human.
if false, the write to the provided file stating that the sample is not human.  

else
set the variable "percent_average" equal to the mean of the percentages contained in percentage_list
write out a brief sentence detailing the percent similarity of the samples  

define the function "human_percentage_blast"
set the variable examplefasta to a given fasta file
set the variable humanfasta to a file containing the human genome
set the variable outputfile to a markdown file that will contain the results of the alignment.  

access the outputfile and designate it as out_handle  

access the file link containing the human genome and designate it as huma.
using gzip, open the human genome file (huma) while reading its contents as text, and save it as the variable "human"  

