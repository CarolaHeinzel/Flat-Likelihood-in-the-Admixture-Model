# Flat-Likelihood-in-the-Admixture-Model

Here, we just describe briefly what the data in this repository is for and how to apply EMALAM. For detailed information, please read the documentation.<br>
# Files in the Repository

This repository contains <br>

EMALAM_incl_J.py: EMALAM, code that the user has to run  <br>
CEU_IBS_TSI_enhanced_corr_f: Example Output of STRUCTURE <br>
Extract_P_J_arbitrary.R: Extract the estimated p in the output of STRUCTURE for J \geq 3 <br>
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
We used python 3.8.3 and <br>

sympy: 1.12 <br>
numpy: 1.23.5 <br>
pandas: 1.4.4 <br>
scipy: 1.10.0 <br>

To install the correct python version, please follow the instrucions at

https://www.python.org/downloads/release/python-383/.  <br>

To install these packages with the correct versions, you can use <br>

python3.8.3 pip install scipy==1.10.0 <br>


and to see whether you where successful, type: <br>
import scipy <br>
print(scipy.__version__) <br>

Alternativly, you can use <br>

conda create -n myenv python=3.8.3 <br>
conda activate myenv <br>
conda install scipy=1.10.0 <br>

to use the correct version of python and to have the correct version of scipy. Analogously, you can install sympy, numpy and pandas. <br>

To run the code, please write <br>

python EMALAM_incl_J.py q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J p_K2_P1 q_K2_P1 "P2" <br>

in the command line. There will be the two output files p_K2_P1_0 q_K2_P1_0 and q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J are the input files. 

# User Choices
Researchers can use different functions to maximize/minimize: <br>
(I) Maximize and minimize the IA for every individuals in every population ("P1") <br>
(II) Maximize the admixture ("P2") <br>
(III) Minimize the admixture ("P3") <br>
(IV) Maximize the ancestry in a specific population ("P4") <br>
(V) Minimize the ancestry in a specific population ("P5") <br>

To choose (I) you can use <br>
python EMALAM_incl_J.py q_CEU_IBS_TSI_K2.txt p_CEU_IBS_TSI_K2 p_CEU_IBS_TSI_K2_J q_K2_P1 p_K2_P1 "P1" <br>

Then, you will get 2K output files for the IAs,  called q_K2_P1_0,  q_K2_P1_1, ..., q_K2_P1_2K and 2K files for the allele frequencies,  called p_K2_P1_0,  p_K2_P1_1, ..., p_K2_P1_2K. <br>

For the other four opportunities, you have to change the last variable of the command line to "P2", "P3", "P4" or "P5". This consequeces two output files, q_K2_P1_0 and  p_K2_P1_0. <br>


You can change the input files to the input files with your own data. Therefore, the same format is required. Specifically, the rows of the input file for the estimated IAs (q_CEU_IBS_TSI_K2.txt) are the individuals and the columns the populations. In p_CEU_IBS_TSI_K2, where the estimated allele frequencies are listed, the rows are the markers and the columns are again the populations. Here, we only need the number of alleles at this marker - 1 rows per marker, e.g. for the output <br>

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



To get the correct format, you can apply the R function Extract_P_J_arbitrary.R (for p_CEU_IBS_TSI_K2_J and p_CEU_IBS_TSI_K2_J) and the R file Extract_p_q.R for q and for cases with only bi-allelic markers. <br>

If you only have bi-allelic markers, type "-1" instead of p_CEU_IBS_TSI_K2_J. 

# Depiction of the Results

To depict the results, you can use the Code Create_Figures.py. This plots the estimated IAs as a bar plot as usually done with STURTCURE results. An example of the output can be found in Documentation_EMALAM.pdf. However, you should apply a software like pong to avoid that your results contain label switching. 








