Import biopython methods  <br>

initialize an empty list, percentage_list.  <br>
  
define the function "add_to_list"<br>
  
if the iteration is not the final iteration, then proceed<br>
add the percentage rounded to 3 decimals to percentage_list.<br>
write this percentage out to a provided file.<br>
checks if the percentage value is greater than or equal to 95%<br>



if true, then write to the provided file stating that the sample is human.<br>
if false, the write to the provided file stating that the sample is not human.<br> 

else<br>
set the variable "percent_average" equal to the mean of the percentages contained in percentage_list<br>
write out a brief sentence detailing the percent similarity of the samples  <br>

define the function "human_percentage_blast"<br>
set the variable examplefasta to a given fasta file<br>
set the variable humanfasta to a file containing the human genome<br>
set the variable outputfile to a markdown file that will contain the results of the alignment.<br>  

access the outputfile and designate it as out_handle  <br>

access the file link containing the human genome and designate it as huma.<br>
using gzip, open the human genome file (huma) while reading its contents as text, and save it as the variable "human"  <br>
Initialize a for loop to parse through the human genome<br>
write the sequence number, id, and sequence record into the provided document<br>

set the variable "aligner" to the value of the method "PairwiseAligner"<br>
set up the aligner settings<br>

Initialize variables to calculate percentage <br>
add percentage to percentage_list <br> 

set a loop limit <br>
