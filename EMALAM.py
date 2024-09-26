import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
import re

import EMALAM_incl_comments_commandline_K2 as em
import Extrahieren_p_q_python as ex
st.set_page_config(page_title="EMALAM", page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

st.header("EMALAM")

cont0 = cont1 = cont2 = False

# Initialize session state variables
if "uploaded_q_file" not in st.session_state:
    st.session_state.uploaded_q_file = None
if "uploaded_p_file" not in st.session_state:
    st.session_state.uploaded_p_file = None
if "poss" not in st.session_state:
    st.session_state.poss = "P4" # Initial Value
if "uploaded_pJ_file" not in st.session_state:
    st.session_state.uploaded_pJ_file = None

def correct_format(data_q):
    data_q = data_q.values.tolist()
    K = len(data_q[0])  # Number of populations
    q_alle = []
    for i in range(K):
        q1_vec = [float(subliste[i]) for subliste in data_q]
        q_alle.append(q1_vec)
    return q_alle
    
    
# Select input data type
data_type = st.selectbox(
    "Select the type of input data:",
    options=["STRUCTURE Output", "Other"],  # Add other options as needed
    index=0  # Set the default selected option to "STRUCTURE"
)

# File uploaders based on input type
if data_type == "STRUCTURE Output":
    st.session_state.uploaded_q_file = st.file_uploader("STRUCTURE output file")
else:
    col1, col2, col3 = st.columns([1,1,1])
    st.session_state.uploaded_q_file = col1.file_uploader("STRUCTURE output file for individual ancestries (q)")
    st.session_state.uploaded_p_file = col2.file_uploader("STRUCTURE output file for allele frequencies (p)")
    st.session_state.uploaded_pJ_file = col3.file_uploader("STRUCTURE output file for allele frequencies (p) for K >= 3")
    markernames = st.checkbox("Markernames", value=False)
    individual_names = st.checkbox("Individual names", value=False)
    if markernames:
    	st.session_state.markernames_file = st.file_uploader("Upload Markernames file")

    if individual_names:
    	st.session_state.individual_names_file = st.file_uploader("Upload Individual Names file")

# Define options for the selectbox
options = [f"P{i+1}" for i in range(5)]

# Function to format the display of the selectbox options
def format_func(poss):
    if poss == "P1":
        return "(I) Maximize individual admixture of one individual"
    elif poss == "P2":
        return "(II) Maximize the average entropy (favours admixed individuals)"
    elif poss == "P3":
        return "(III) Minimize the average entropy (favours non-admixed individuals)"
    elif poss == "P4":
        return "(IV) Maximize individual admixture of one population"
    elif poss == "P5":
        return "(IV) Minimize individual admixture of one population"

# Selectbox for choosing the function, linked to session state 'poss'
poss = st.selectbox(
    "Function which is optimized", 
    options, 
    index=None, 
    format_func=format_func, 
    key="poss",  # Bind to session state key 'poss'
    help="See Documentation", 
    placeholder="Choose an option", 
    disabled=False, 
    label_visibility="visible"
)

if(data_type != "STRUCTURE Output"):
# Check if files are uploaded before processing them
	if st.session_state.uploaded_p_file is not None:
    		data_p = pd.read_csv(st.session_state.uploaded_p_file, delimiter=" ", header=None)
    		M = data_p.shape[0]    
    		K = data_p.shape[1]

	if st.session_state.uploaded_pJ_file is not None:
    		data_pJ = pd.read_csv(st.session_state.uploaded_pJ_file, delimiter=" ", header=None)
    		M = data_pJ.shape[0]    
    		K = data_pJ.shape[1]

	if st.session_state.uploaded_q_file is not None:
    		data_q = pd.read_csv(st.session_state.uploaded_q_file, delimiter=" ", header=None)
    		N = data_q.shape[5]
    		st.write(data_q)
    		st.write(data_p)
    		st.write(data_pJ)

else:
	K = st.number_input("Number of Populations:", min_value=1, step=1)

	if st.session_state.uploaded_q_file is not None:
	
	  	#with open(file_path, 'r') as file:
	  	uploaded_file = st.session_state.uploaded_q_file.readlines()

	  	data_pJ, data_p, marker_names = ex.read_table_data(uploaded_file, K)
	  	data_q, individual_names = ex.extract_q(uploaded_file, K)
	  	M = data_p.shape[0]  
	  	N = data_q.shape[0]
	  	st.write(data_q, data_p, data_pJ)

selected_individual = 0   
if poss == "P1":
    individual_options = range(0, N) 
    # To do: in algo_final noch einbauen
    selected_individual = st.selectbox("Which individual is maximized?", individual_options)
if poss == "P4" or poss == "P5":
    individual_options = range(0, K) 
    k_specific = st.selectbox("Which population should be considered?", individual_options)
submit = st.button("Submit")
if submit:
    # Entfernen nicht-polymorpher Loci
    data_p = data_p[~(data_p == 1).all(axis=1)]
    pJ = data_pJ
    simi = 0
    k_specific = 0
    n_trials = 10  # Anzahl der verschiedenen Anfangswerte für die Minimierungsfunktion

    # Überprüfen der Form der Daten
    if data_q.shape[1] != K:
        st.write("Die Anzahl der Spalten in der q-Datei muss der Anzahl der Spalten in der p-Datei (Anzahl der Populationen) entsprechen.")
        st.rerun()

    data_p = pd.DataFrame(data_p)
    data_q_out, data_p_out = em.algo_final(data_q, data_p, K, poss, simi, k_specific, data_pJ, n_trials, selected_individual)

    # Eingabedaten anzeigen
    with st.expander('Input Data'):
        cols = st.columns([K] + [1 for _ in range(K)])
        with cols[0]:
            st.subheader("Daten q")
            st.bar_chart(data_q, horizontal=True)
        with cols[1]:
            st.subheader("Daten p")
            st.bar_chart(data_p, horizontal=True)

    # Individuelle Admixture (q) darstellen
    with st.expander('Individuelle Admixture (q)'):
        cols = st.columns(len(data_q_out))
        print(data_q_out)
        for i, d in enumerate(data_q_out):
            with cols[i]:
                if K == 2:  # Überprüfen, ob K gleich 2 ist
                    st.write(d["data"])
                    data_q = d["data"]
                    colors = plt.cm.tab10(np.linspace(0, 1, data_q.shape[1]))

                    q_alle = correct_format(data_q)  # Korrekte Formatierung der Daten
                    q_alle_1 = np.transpose(q_alle) 

                    # Erstellen des gestapelten Balkendiagramms
                    fig, ax = plt.subplots(figsize=(20, 8))  # Breite und Höhe nach Bedarf anpassen
                    bottom = np.zeros(q_alle_1.shape[0])  # Basislinie für das gestapelte Diagramm

                    for j in range(q_alle_1.shape[1]):  # Über jede Population iterieren
                        ax.bar(range(q_alle_1.shape[0]), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                        bottom += q_alle_1[:, j]  # Aktualisieren der Basislinie

                    individual_names = [f'Individuum {k + 1}' for k in range(q_alle_1.shape[0])]  # Individuelle Namen
                    ax.set_xticks(range(len(individual_names)))  # Setze die Positionen der x-Achsen-Ticks
                    ax.set_xticklabels(individual_names, rotation=90, ha='right')  # Setze die benutzerdefinierten Labels und drehe sie

                    ax.set_xlabel('Individuals')
                    ax.set_ylabel('Estimated IAs')
                    ax.set_ylim(0, 1)
                    ax.legend(title='Populationen')  # Legende hinzufügen

                    # Grafik in Streamlit anzeigen
                    st.pyplot(fig)
                else:  # Überprüfen, ob K gleich 2 ist
                
                    st.write(data_q_out[0])
                    data_q = data_q_out[0]
                    colors = plt.cm.tab10(np.linspace(0, 1, K))

                    q_alle = data_q
                    q_alle_1 = np.transpose(q_alle) 
                    print(q_alle_1)
		     
                    # Erstellen des gestapelten Balkendiagramms
                    fig, ax = plt.subplots(figsize=(20, 8))  # Breite und Höhe nach Bedarf anpassen
                    bottom = np.zeros(N)  # Basislinie für das gestapelte Diagramm

                    for j in range(K):  # Über jede Population iterieren
                        ax.bar(range(N), q_alle_1[:, j], bottom=bottom, color=colors[j], label=f'Population {j + 1}')
                        print(j, len(q_alle_1[:,j]), len(bottom))
                        for l in range(len(bottom)):
                        	bottom[l] += q_alle_1[l, j]  

                    individual_names = [f'Individuum {k + 1}' for k in range(q_alle_1.shape[0])]  # Individuelle Namen
                    ax.set_xticks(range(len(individual_names)))  # Setze die Positionen der x-Achsen-Ticks
                    ax.set_xticklabels(individual_names, rotation=90, ha='right')  # Setze die benutzerdefinierten Labels und drehe sie

                    ax.set_xlabel('Individuals')
                    ax.set_ylabel('Estimated IAs')
                    ax.set_ylim(0, 1)
                    ax.legend(title='Populationen')  # Legende hinzufügen

                    # Grafik in Streamlit anzeigen
                    st.pyplot(fig)

    with st.expander(f'Graphical representation of p (allele frequencies)'):
        cols = st.columns([1 for i in data_p_out])        
        i = 0
        for d in data_p_out:
            with cols[i]:
                # st.write(d["data"])
                # st.write(d["extension"])
                co = st.columns([1 for i in range(K)])
                for j in range(K):
                    with co[j]:
                        #p_df_loc = pd.DataFrame({'0' : d["data"][j], '1': 1-d["data"][j]})
                        #st.bar_chart(p_df_loc, horizontal=True)
                        #st.bar_chart(d["data"][j], horizontal=True)
                        st.write(data_p_out)
                i = i+1






