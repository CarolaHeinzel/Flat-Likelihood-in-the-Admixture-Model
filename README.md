# Usage

After installing dependencies using `pip install -r requirements.txt`, either use the command line tool, e.g.
```
python emalam.py --structure_filename Example_Input/CEU_IBS_TSI_enhanced_corr_K3_f --out out.txt out.json --fun entropy --min 
```
or the graphical version
```
streamlit run EMALAM.py
```
which is displayed in your browser. Another online version is available under https://flat-likelihood-in-the-admixture-model-p4vvpwxjdtteu3prahatyj.streamlit.app/.


# Flat-Likelihood-in-the-Admixture-Model

We briefly describe how to apply emalam. For detailed information, please read the documentation.<br>
# Files in the Repository

This repository contains <br>

EMALAM_incl_comments_commandline_K2.py: EMALAM, code that the user has to run  <br>
CEU_IBS_TSI_enhanced_corr_f: Example Output of STRUCTURE <br>
Extract_P_J_arbitrary.R: Extract the estimated p in the output of STRUCTURE for J > 2 <br>
Extract_p_q.R: Extract the estimated p and estimated q in the output of STRUCTURE for J = 2 <br>
Documentation_EMALAM.pdf: Explains how to use EMALAM. <br>
Create_Plot.py: Depicts the estimated IAs <br>

Additionally, the data 

p_CEU_IBS_TSI_K2 <br>
p_CEU_IBS_TSI_K2_J <br>
q_CEU_IBS_TSI_K2.txt <br>

are the three input files (as an example) and

p_K2_P1_0.txt <br>
q_K2_P1_0.txt <br>

are an example output for possibility I. The definition of possibility I is described below.

# Run the Code
We used python 3.8.3 and the packages as mentioned in requirements.txt. To install them, type `pip install -r requirements.txt`. We have tested the code under python 3.8 and 3.12. 
 <br>
 
To run the code, please write <br>

python EMALAM_incl_comments_commandline_K2.py q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J p_K2_P1 q_K2_P1 "P2" <br>

in the command line. There will be the two output files p_K2_P1_0, q_K2_P1_0 and q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J are the input files. We describe the meaning of "P2" below.

# User Choices
Researchers can use different functions to maximize/minimize: <br>
(I) Maximize and minimize the estimated IA for every individual in every population ("P1") <br>
(II) Maximize the estimated admixture ("P2") <br>
(III) Minimize the estimated admixture ("P3") <br>
(IV) Maximize the estimated ancestry in a specific population ("P4") <br>
(V) Minimize the estimated ancestry in a specific population ("P5") <br>

To choose (I) you can use <br>
python EMALAM_incl_comments_commandline_K2.py q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J q_K2_P1 p_K2_P1 "P1" <br>

Then, you will get 2K output files for the estimated IAs,  called q_K2_P1_0.txt,  q_K2_P1_1.txt, ..., q_K2_P1_2K.txt and 2K files for the estimated allele frequencies,  called p_K2_P1_0.txt,  p_K2_P1_1.txt, ..., p_K2_P1_2K.txt. <br>

For the other four opportunities, you have to change the last variable of the command line to "P2", "P3", "P4" or "P5". This consequeces two output files, q_K2_P1_0.txt and  p_K2_P1_0.txt. If you chose "P4" or "P5", an additional input, i..e. the population that should be considered, is required. To consider population i = 0, 1, ..., K, you type <br>

python EMALAM_incl_comments_commandline_K2.py q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J q_K2_P1 p_K2_P1 "P4" --k_specific i <br>



You can change the input files to the input files with your own data. Therefore, the same format is required. Specifically, the rows of the input file for the estimated IAs (q_CEU_IBS_TSI_K2.txt) are the individuals and the columns the populations. In p_CEU_IBS_TSI_K2, where the estimated allele frequencies are listed, the rows are the markers and the columns are again the populations. Here, we only need (the number of alleles at this marker - 1) rows per marker, e.g. for the output <br>

Locus 1 : <br>
2 alleles <br>
0.0% missing data <br>
   1   (0.391) 0.391 0.391 0.269 <br>
   0   (0.609) 0.609 0.609 0.731 <br>

  we need the line <br>

0.609 0.609 0.731  <br>

or <br>

0.391 0.391 0.269 <br>

as input file. <br>

The rows of p_CEU_IBS_TSI_K2_J also represent the markers and the columns represent the populations. However, here we sum up the allele frequencies of every marker with more than two alleles. For example, the output of STRUCTURE might be <br>

Locus 1 :  <br>
2 alleles <br>
0.0% missing data <br>
   1   (0.391) 0.391 0.391 0.269  <br>
   0   (0.609) 0.609 0.609 0.731  <br>

Locus 2 :  <br>
2 alleles <br>
0.0% missing data <br>
   0   (0.240) 0.640 0.840 0.576  <br>
   1   (0.160) 0.160 0.160 0.124   <br>
   2   (0.600) 0.200 0.000 0.300  <br>
Then, the file p_CEU_IBS_TSI_K2_J contains

0.360 0.160 0.424   <br>

To get the correct format, you can apply the R function Extract_P_J_arbitrary.R (for p_CEU_IBS_TSI_K2_J and p_CEU_IBS_TSI_K) and the R file Extract_p_q.R for q and for cases with only bi-allelic markers. <br>

If you only have bi-allelic markers, type "-1" instead of p_CEU_IBS_TSI_K2_J. 

# Depiction of the Results

To depict the results, you can use the Code Create_Figures.py. This plots the estimated IAs as a bar plot as usually done with STURTCURE results. An example of the output can be found in Documentation_EMALAM.pdf. However, you should apply a software like pong to avoid that your results contain label switching. 








