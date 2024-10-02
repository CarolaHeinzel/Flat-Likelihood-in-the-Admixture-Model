import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
import re

import EMALAM_incl_comments_commandline_K2 as em
import Extrahieren_p_q_python as ex

st.set_page_config(page_title="EMALAM", page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

# Initialize session state variables
if "uploaded_structure_file" not in st.session_state:
    st.session_state.uploaded_structure_file = None
if "uploaded_q_file" not in st.session_state:
    st.session_state.uploaded_q_file = None
if "uploaded_p_file" not in st.session_state:
    st.session_state.uploaded_p_file = None
if "poss" not in st.session_state:
    st.session_state.poss = "P4" # Initial Value
if "uploaded_pJ_file" not in st.session_state:
    st.session_state.uploaded_pJ_file = None
if "data_uploaded" not in st.session_state:
    st.session_state.data_uploaded = None

def load_default_file():
    default_file_path = 'Example_Output_STRUCTURE'  
    with open(default_file_path, 'rb') as f:
        return f.readlines()

default_file_path_q = 'q_CEU_IBS_TSI_K2.txt'
default_file_path_p = 'p_CEU_IBS_TSI_K2' 
default_file_pJ = 0 

def load_default_file_notSTRUCTURE(default_file_path_q, default_file_path_p, default_file_pJ):
    data_p = pd.read_csv(default_file_path_p, delimiter=" ", header=None)
    data_q = pd.read_csv(default_file_path_q, delimiter=" ", header=None)
    if default_file_pJ:
        data_J = pd.read_csv(default_file_pJ, delimiter=" ", header=None)
    else:
         data_J = 0     
    return data_q, data_p, data_J

def correct_format(data_q):
    data_q = data_q.values.tolist()
    K = len(data_q[0])  # Number of populations
    q_alle = []
    for i in range(K):
        q1_vec = [float(subliste[i]) for subliste in data_q]
        q_alle.append(q1_vec)
    return q_alle
    
# Determine K for the STRUCTURE File    
def extract_population_count(file1):
    output_text = ''.join([line.decode('utf-8') for line in file1])
    match = re.search(r"(\d+)\s+populations assumed", output_text)

    if match:
        return int(match.group(1))  
    else:
        return None 

st.header("EMALAM")
st.write("See xxx for a reference what this webpage is about.")
st.write("You can use an output file of STRUCTURE, which will then be sorted into $q$- and $p$-matrices.")
st.write("For the other case, you must upload two files (for bi-allelic loci) or three files (of some alleles are multi-allelic).:\n  * The $q$-file of individual admixtures;\n  * The $p$-file of allele frequencies in all populations;\n  * Some J file??")

default_options = ["Example [STRUCTURE output data](github.com)", "Example [other data](github.com)", "Upload own data"]
default = st.radio(
    "Which file(s) should be used?",
    options = default_options,
    index = 0,
)

if default == default_options[0]:
    use_structure_file = True
    ## read default files
elif default == default_options[1]:
    use_other_file = True
    ## read default files    
elif default == default_options[2]:
    # Select input data type
    data_type_options = ["STRUCTURE Output", "Other"]
    data_type = st.selectbox(
        "Select the type of input data:",
        options=data_type_options,  # Add other options as needed
        index=0  # Set the default selected option to "STRUCTURE"
    )
    use_structure_file = (data_type == data_type_options[0])
    use_other_file = (data_type == data_type_options[1])
    
    # File uploaders based on input type
    if data_type == data_type_options[0]:
        use_structure_file = True
        st.session_state.uploaded_structure_file = st.file_uploader("STRUCTURE output file")
    else:
        use_other_file = True
        col1, col2, col3 = st.columns([1,1,1])
        st.session_state.uploaded_q_file = col1.file_uploader("STRUCTURE output file for individual ancestries (q)")
        st.session_state.uploaded_p_file = col2.file_uploader("STRUCTURE output file for allele frequencies (p)")
        st.session_state.uploaded_pJ_file = col3.file_uploader("STRUCTURE output file for allele frequencies (p) for J >= 3")



    #markernames = st.checkbox("Markernames", value=False)
    #individual_names = st.checkbox("Individual names", value=False)
    #if markernames:
    #	st.session_state.markernames_file = st.file_uploader("Upload Markernames file")

    #if individual_names:
    #	st.session_state.individual_names_file = st.file_uploader("Upload Individual Names file")



if st.session_state.data_uploaded == True:
    # Selectbox for choosing the function, linked to session state 'case'
    case_options = ["(I) Minimize/Maximize average contribution from population i", "(II) Minimize/Maximize average entropy"]
    # Add other options as needed

    case = st.selectbox(
        "Function which is optimized:",
        options=case_options, 
        index=0,  # Set the default selected option to "entropy"
        help = "Maximal entropy favours admixed individuals, minimal entropy favours non-admixed individuals."
    )


if data_type != data_type_options[0]:
    if st.session_state.uploaded_p_file is not None:
        data_p = pd.read_csv(st.session_state.uploaded_p_file, delimiter=" ", header=None)
        M = data_p.shape[0]
        K = data_p.shape[1]
        m

    if st.session_state.uploaded_pJ_file is not None:
        data_pJ = pd.read_csv(st.session_state.uploaded_pJ_file, delimiter=" ", header=None)
        K = data_pJ.shape[1]

    if st.session_state.uploaded_q_file is not None:
        data_q = pd.read_csv(st.session_state.uploaded_q_file, delimiter=" ", header=None)
        N = data_q.shape[0]
        individual_names = np.linspace(1, N, N)
        st.write(data_q)
        st.write(data_p)
        st.write(data_pJ)
    else:
        data_q, data_p, data_pJ = load_default_file_notSTRUCTURE()
        M = data_p.shape[0]
        K = data_p.shape[1]
        N = data_q.shape[0]
        st.write(data_q)
        st.write(data_p)
        st.write(data_pJ)
        individual_names = np.linspace(1, N, N) 
        marker_names = np.linspace(1, M, M)

    		
else:

	if st.session_state.uploaded_q_file is not None:
	
	  	uploaded_file = st.session_state.uploaded_q_file.readlines()
	  	K = extract_population_count(uploaded_file)

	  	data_pJ, data_p, marker_names, p_all = ex.read_table_data(uploaded_file, K)
	  	data_q, individual_names = ex.extract_q(uploaded_file, K)
	  	p_all = ex.(uploaded_file, K)
	  	print(p_all)
	  	M = data_p.shape[0]  
	  	N = data_q.shape[0]
	  	st.write(data_q, data_p, data_pJ)
	  	print(type(data_pJ))
	  	if(type(data_pJ) != np.ndarray):
	  		if(data_pJ == None):
	  			data_pJ = 0
	else:
   	 	uploaded_file = load_default_file()
   	 	K = extract_population_count(uploaded_file)
   	 	data_pJ, data_p, marker_names, p_all = ex.read_table_data(uploaded_file, K)
   	 	print(p_all)
   	 	data_q, individual_names = ex.extract_q(uploaded_file, K)

   	 	M= data_p.shape[0] 
   	 	N = data_q.shape[0]
   	 	st.write(data_q, data_p, data_pJ)
   	 	if(type(data_pJ) != np.ndarray):
   	 		if(data_pJ == None):
   	 			data_pJ = 0

selected_individual = 0   
if poss == "P1":
    individual_options = range(0, N) 
    selected_individual = st.selectbox("Which individual should be considered?", individual_options)
if poss == "P4" or poss == "P5":
    individual_options = range(0, K) 
    k_specific = st.selectbox("Which population should be considered?", individual_options)
submit = st.button("Submit")

if submit:
    with st.spinner('The computer is calculating...'):


	    data_p = data_p[~(data_p == 1).all(axis=1)]
	    pJ = data_pJ
	    simi = 0
	    k_specific = 0
	    n_trials = 10



	    data_p = pd.DataFrame(data_p)
	    data_q_out, data_p_out = em.algo_final(data_q, data_p, K, poss, simi, k_specific, data_pJ, n_trials, selected_individual, data_type)

    with st.expander('Input Data'):
        cols = st.columns([K] + [1 for _ in range(K)])
        with cols[0]:
            st.subheader("Data q")
            data_q = data_q
            colors = plt.cm.tab10(np.linspace(0, 1, data_q.shape[1]))

            q_alle = correct_format(data_q)
            q_alle_1 = np.transpose(q_alle)

            fig, ax = plt.subplots(figsize=(20, 8))
            bottom = np.zeros(q_alle_1.shape[0])

            for j in range(q_alle_1.shape[1]):
                ax.bar(range(q_alle_1.shape[0]), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                bottom += q_alle_1[:, j]

            ax.set_xticks(range(len(individual_names)))
            #print(individual_names)

            ax.set_xticklabels(individual_names, rotation=90, ha='right')
            ax.set_xlabel('Individuals')
            ax.set_ylabel('Estimated IAs')
            ax.set_ylim(0, 1)
            ax.legend(title='Populations')

            st.pyplot(fig)
        with cols[1]:
            st.subheader("Data p")
            data_q = data_p
            colors = plt.cm.tab10(np.linspace(0, 1, data_q.shape[1]))

            q_alle = correct_format(data_q)
            q_alle_1 = np.transpose(q_alle)

            fig, ax = plt.subplots(figsize=(20, 8))
            bottom = np.zeros(q_alle_1.shape[0])

            for j in range(q_alle_1.shape[1]):
                ax.bar(range(q_alle_1.shape[0]), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                bottom += q_alle_1[:, j]

            ax.set_xticks(range(len(marker_names)))
            ax.set_xlabel('Markers')
            ax.set_ylabel('Estimated Allele Frequencies')
            ax.set_ylim(0, K)
            ax.legend(title='Populations')

            st.pyplot(fig)

    # Individuelle Admixture (q) darstellen
    with st.expander('Individual Admixture (q)'):
        cols = st.columns(len(data_q_out))

        for i, d in enumerate(data_q_out):
            st.write(data_q_out)
        
            with cols[i]:
                if K == 2 and poss != "P2" and poss != "P3":  
                    #st.write(d["data"])
                    data_q = d["data"]
                    colors = plt.cm.tab10(np.linspace(0, 1, data_q.shape[1]))

                    q_alle = correct_format(data_q)  # Korrekte Formatierung der Daten
                    q_alle_1 = np.transpose(q_alle) 

                    fig, ax = plt.subplots(figsize=(20, 8))  # Breite und Höhe nach Bedarf anpassen
                    bottom = np.zeros(q_alle_1.shape[0])  # Basislinie für das gestapelte Diagramm

                    for j in range(q_alle_1.shape[1]):  # Über jede Population iterieren
                        ax.bar(range(q_alle_1.shape[0]), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                        bottom += q_alle_1[:, j]  # Aktualisieren der Basislinie

                    ax.set_xticks(range(len(individual_names)))  # Setze die Positionen der x-Achsen-Ticks
                    ax.set_xticklabels(individual_names, rotation=90, ha='right')  # Setze die benutzerdefinierten Labels und drehe sie

                    ax.set_xlabel('Individuals')
                    ax.set_ylabel('Estimated IAs')
                    ax.set_ylim(0, 1)
                    ax.legend(title='Populationen')  

                    st.pyplot(fig)
                else:  
                
                    st.write(data_q_out[0])
                    data_q = data_q_out[0]
                    colors = plt.cm.tab10(np.linspace(0, 1, K))

                    q_alle = data_q
                    q_alle_1 = np.transpose(q_alle) 
		     
                    fig, ax = plt.subplots(figsize=(20, 8))  
                    bottom = np.zeros(N)  

                    for j in range(K):  
                        ax.bar(range(N), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                        for l in range(len(bottom)):
                        	bottom[l] += q_alle_1[l, j]  
                    ax.set_xticks(range(len(individual_names)))  
                    ax.set_xticklabels(individual_names, rotation=90, ha='right') 

                    ax.set_xlabel('Individuals')
                    ax.set_ylabel('Estimated IAs')
                    ax.set_ylim(0, 1)
                    ax.legend(title='Populations')  

                    st.pyplot(fig)

    with st.expander(f'Allele Frequencies (p) '):
        cols = st.columns([1 for i in data_p_out])        
        i = 0
        st.write(data_p_out)

        for d in data_p_out:
            with cols[i]:

                if(K == 2 and poss != "P2" and poss != "P3"):
                

                    data_q = d["data"]
                    colors = plt.cm.tab10(np.linspace(0, 1, data_q.shape[1]))

                    q_alle = correct_format(data_q)  
                    q_alle_1 = np.transpose(q_alle) 

                    fig, ax = plt.subplots(figsize=(20, 8))  
                    bottom = np.zeros(q_alle_1.shape[0])  
		    
                    for j in range(q_alle_1.shape[1]): 
                        ax.bar(range(q_alle_1.shape[0]), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                        bottom += q_alle_1[:, j]  

                    #ax.set_xticks(range(len(marker_names)))  
                    #ax.set_xticklabels(marker_names, rotation=90, ha='right')  

                    ax.set_xlabel('Markers')
                    ax.set_ylabel('Estimated Allele Frequencies')
                    ax.set_ylim(0, K)
                    ax.legend(title='Populations')  

                    st.pyplot(fig)                
                else:  
                
                    st.write(data_p_out[0])
                    data_q = data_p_out[0]
                    colors = plt.cm.tab10(np.linspace(0, 1, K))

                    q_alle = data_q
                    q_alle_1 = np.transpose(q_alle) 
		     
                    fig, ax = plt.subplots(figsize=(20, 8))  
                    bottom = np.zeros(M)  
                    bar_width = 0.15

                    for j in range(K):  
                    	x_positions = np.arange(M) + j * bar_width  # Position für Population j
                    	ax.bar(x_positions, q_alle_1[:, j], width=bar_width, color=colors[j], label=f'Population {j + 1}')
			    
                    ax.set_xticks([])  

                    ax.set_xlabel('Markers')
                    ax.set_ylabel('Estimated Allele Frequencies')
                    ax.set_ylim(0, 1)
                    ax.legend(title='Populations')  

                    st.pyplot(fig)



                i = i+1
    if (data_type == "STRUCTURE Output"):

    	with st.expander(f'Marker Names and Individual Names '):

                st.write(marker_names, individual_names)


