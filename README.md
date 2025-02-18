# Usage

After installing dependencies using `pip install -r requirements.txt`, either use the command line tool, e.g.
```
python emalam.py --structure_filename Example_Input/CEU_IBS_TSI_enhanced_corr_K3_f --out out.txt out.json --fun entropy --min 
```
or (for ADMIXTURE output)
```
python emalam.py --hatq_filename combined_filtered.3.Q --hatp_filename combined_filtered.3.P --out outEURmax.Q outEURmax.P --fun size --pop 2 --max
```
or the graphical version
```
streamlit run EMALAM.py
```

which is displayed in your browser. Another online version is available under [streamlit](https://flat-likelihood-in-the-admixture-model-p4vvpwxjdtteu3prahatyj.streamlit.app/).
The code is tested with python 3.12.
 
# User Choices

The choices of the user are described in EMALAM\_Manual.pdf.

# Depiction of the Results

To depict the results, you can use the Code Create_Figures.py. This plots the estimated IAs as a bar plot as usually done with STURTCURE results. 








