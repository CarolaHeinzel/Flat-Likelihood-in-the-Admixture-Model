# Flat-Likelihood-in-the-Admixture-Model

Here, we just describe briefly what the data in this repository is for. For detailed information, please read the documentation.<br>

Researchers can use different functions to maximize/minimize: <br>
(I) Maximize and minimize the IA for every individuals in every population <br>
(II) Maximize the admixture<br>
(III) Minimize the admixture <br>
(IV) Maximize the ancestry in a specific population <br>
(V) Minimize the ancestry in a specific population <br>


This repository contains <br>

EMALAM_incl_J.py: EMALAM, code that the user has to run  <br>
CEU_IBS_TSI_enhanced_corr_f: Example Output of STRUCTURE <br>
Extract_P_J_arbitrary.R: Extract the estimated p in the output of STRUCTURE for J \geq 3 <br>
Extract_p_q.R: Extract the estimated p and estimated q in the output of STRUCTURE for J = 2 <br>
Documentation_EMALAM: Explains how to use EMALAM. <br>
Create_Plot.py: Depicts the estimated IAs <br>

Additionally, the data 

p_CEU_IBS_TSI_K2 <br>
p_CEU_IBS_TSI_K2_J <br>
q_CEU_IBS_TSI_K2.txt <br>

are the three input files (as an example) and

p_K2_P1_0.txt <br>
q_K2_P1_0.txt <br>

are an example output for possibility I.

We used python 3.8.3 and <br>

sympy: 1.12 <br>
numpy: 1.23.5 <br>
pandas: 1.4.4 <br>
scipy: 1.10.0 <br>

To install these packages, you can use <br>

python3.8.3 pip install scipy==1.10.0 <br>

and to see whether you where successful, type:
import scipy
print(scipy.__version__) <br>

Alternativly, you can use <br>

conda create -n myenv python=3.8.3 <br>
conda activate myenv <br>
conda install scipy=1.10.0 <br>

to use the correct version of python and to have the correct version of scipy. Analogously, you can install sympy, numpy and pandas.

