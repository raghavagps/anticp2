# Anticp2: Prediction, Design and scan of anticancer prptides
AntiCP 2.0 is an updated version of AntiCP, developed to predict and design anticancer peptides with high accuracy. This study utilize largest possible dataset of anticancer and non-anticancer peptides. Main dataset consists of experimentally validated 861 anticancer peptides and 861 non-anticancer or validated antimicrobial peptides. Alternate dataset comprises of 970 anti-cancer peptides and 970 non-anticancer peptides (randomly pickup from Swiss-Prot). 
# Reference
Agrawal P., Bhagat D., Mahalwal M., Sharma N., Raghava G.P.S. (2020), AntiCP 2.0: an updated model for predicting anticancer peptides, <a href="https://doi.org/10.1093/bib/bbaa153">Briefings in Bioinformatics, bbaa153</a> 
# Web Server
https://webs.iiitd.edu.in/raghava/anticp2/

## PIP Installation
PIP version is also available for easy installation and usage of this tool. The following command is required to install the package 
```
pip install anticp2
```
To know about the available option for the pip package, type the following command:
```
anticp2 -h
```
# Installation
git clone https://github.com/raghavagps/anticp2

change dir to anticp2
# Requirement 
It is developed for python3 and require following libraries, these libraries (pandas, numpy, sklearn, pickle-mixin) can be install using following commands
```
pip install pandas    
pip install numpy   
pip install sklearn     
pip install pickle-mixin
```
# Introduction
AntiCP2 is developed for predicting, desiging and scanning antcancer peptides. More information on AntiCP2 is abvailble from its web server http://webs.iiitd.edu.in/raghava/anticp2/ . This page provide information about stnadalone version of AntiCP2. Please read/cite following paper for complete information including algorithm behind AntiCP2.

Agrawal P., Bhagat D., Mahalwal M., Sharma N., and Raghava GPS (2020) AntiCP 2.0: an updated model for predicting anticancer peptides. Briefing in Bioinformatics doi: 10.1093/bib/bbaa153

**Models:** In this program, two models have beeen incorporated for predicting anticancer peptides. Model1 is trained on Anti-Cancer and Anti-Microbial peptides, it is default model. Model2 is trained on Anti-Cancer and Non-Anticalcer (or random peptides) peptides.

**Modules/Jobs:** This program implement three modules (job types); i) Predict: for predictin anticancer peptides, ii) Design: for generating all possible peptides and computing Anti-Cancer potential (score) of peptides, iii) Scan: for creating all possible overlapping peptides of given length (window) and computing Anti-Cancer potential (score) of these overlapping peptides.

**Minimum USAGE:** Minimum ussage is "anticp2 -i peptide.fa" where peptide.fa is a input fasta file. This will predict Anti-Cancer potential of sequence  in fasta format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

**Full Usage:** Following is complete list of all options, you may get these options by "anticp2 -h" 

```
anticp.py [-h] -i INPUT [-o OUTPUT] [-j {1,2,3}] [-t THRESHOLD]
                  [-m {1,2}]
                  [-w {5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29}]
                  [-d {1,2}]
```

**optional arguments:**
```
  -h, --help            show this help message and exit

  -i INPUT, --input INPUT
                        Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code

-o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv

-j {1,2,3}, --job {1,2,3}
                        Job Type: 1:predict, 2:design and 3:scan, by default 1

-t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.5

-m {1,2}, --model {1,2}
                        Model: 1: ACP/AMP, 2: ACP/non-ACP, by default 1

-w {5,6,7,..,30}, --winleng
                        Window Length: 5 to 30 (scan mode only), by default 10

-d {1,2}, --display {1,2}
                        Display: 1:Anticancer peptide, 2: All peptides, by default 1
```

**Input File:** It allow users to provide input in two format; i) FASTA format (standard) and ii) Simple Format. In case of simple format, file should have one one peptide sequence in a single line in single letter code (eg. peptide.seq). Please note in case of predict and design module (job) length of peptide should be upto 50 amino acids, if more than 50, program will take first 50 residues. In case of of scan module, minimum length of protein/peptide sequence should be more than equal to window length (pattern), see peptide.fa . Please note program will ignore peptides having length less than 5 residues (e.g., protein.fa).

**Output File:** Program will save result in CSV format, in case user do not provide output file name, it will be stored in outfile.csv.

**Threshold:** User should provide threshold between 0 and 1, please note score is propotional to anti-cancer potential of peptide.


## AntiCP2 Packakage Files

It contantain following files, brief descript of these files given below
```
INSTALLATION  	: Installations instructions

LICENSE       	: License information

README.md     	: This file provide information about this package

aac_extra_model : Model file required for running Model 2

anticp2.py 	: Main python program 

dpc_extra_model : Model file required for running Model 1 (Default) 

outfile.csv	: Example output file in csv format

peptide.fa	: Example file contain peptide sequenaces in FASTA format

peptide.seq	: Example file contain peptide sequenaces in simple format

protein.fa	: Example file contain protein sequenaces in FASTA format 

```
# Address for contact
In case of any query please contact 
```
Prof. G. P. S. Raghava, Head Department of Computational Biology,            
Indraprastha Institute of Information Technology (IIIT), 
Okhla Phase III, New Delhi 110020 ; Phone:+91-11-26907444; 
Email: raghava@iiitd.ac.in  Web: http://webs.iiitd.edu.in/raghava/
```



