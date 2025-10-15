Import biopython methods  <br>
Import numpy, urllib, and gzip

initialize an empty list, percentage_list.  <br>
  
define a function to add percentages to a list<br>
  
if the iteration is not the final iteration, then proceed<br>
add the percentage rounded to 3 decimals to percentage_list.<br>
write this percentage out to a provided file.<br>
checks if the percentage value is greater than or equal to 95%<br>


if true, then write to the provided file stating that the sample is human.<br>
if false, the write to the provided file stating that the sample is not human.<br> 

else<br>
set the variable "percent_average" equal to the mean of the percentages contained in percentage_list<br>
write out a brief sentence detailing the percent similarity of the samples  <br>

define the function to read the input file<br>
  open and read the input file.<br>

define the function to read the file containing the human genome <br>
first, open the link containing the human genome file <br>
open the file containing the human genome <br>
parse through the data in the fasta format. <br>

define the function for alignment calculations <br>
set up the aligner using the PairwiseAligner <br>
set up the aligner settings<br>

define the main method <br>
have a variable to take an input fasta file <br>
have a variable that contains the human genome <br>
set an output file to write the results into <br>
open the output file with writing permissions <br>
iterate through the human genome <br>
write into the output based on what is read in the human genome <br>
use the alignment calculator function to calculate the percentage <br>
add the percentage to the percentage list, and writes it into the output file <br>
set a read limit so that the program doesn't iterate through the entire human genome <br>

