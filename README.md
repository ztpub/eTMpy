# eTMplant
Endogenous target mimic (eTMs) screen package


Based on the following rigorous rules:

(1) three unpaired nucleotides (bulges) are only permitted at the positions corresponding to the 9th to 12th nucleotides counting from the 5¡¯ end of a miRNA sequence; 

(2) the 2nd to 8th positions at the 5¡¯ of a miRNA sequences must be perfectly pairing; 

(3) the total number of mismatches and G/U pairs within the eTM and miRNA pairing region (excluding the bulge region) should be no more than three and the consecutive mismatches should not exceed two. 


######################################################################

Demo:
 
   run python eTMplant.py Arabidopsis_thaliana.fasta ath_miRNA.txt

Results will be in eTM_py_list.txt

######################################################################


* Version update 0.2
* Add the screen of postive and negative chains together
