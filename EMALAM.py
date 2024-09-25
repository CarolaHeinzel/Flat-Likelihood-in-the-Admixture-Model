import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

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
        N = data_q.shape[0]
        st.write(data_q)
        st.write(data_p)
        st.write(data_pJ)

else:
    K = st.number_input("Number of Populations:", min_value=1, step=1)
    if st.session_state.uploaded_q_file is not None:
	  	#with open(file_path, 'r') as file:
        uploaded_file = st.session_state.uploaded_q_file.readlines()

        data_pJ, data_p = ex.read_table_data(uploaded_file, K)
        data_q = ex.extract_q(uploaded_file, K)
        M = data_p.shape[0]  
        N = data_q.shape[0]

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
    # Remove non-polymorphic loci
    data_p = data_p[~(data_p == 1).all(axis=1)]

    pJ = data_pJ
    simi = 0
    k_specific = 0
    n_trials = 10 # Number of different initial values for the minimization function

    if data_q.shape[1] != K:
        st.write("Number of columns in q-file must equal number of columns in p-file (number of populations)")
        st.rerun()

    data_p = pd.DataFrame(data_p)
    data_q_out, data_p_out  = em.algo_final(data_q, data_p, K, poss, simi, k_specific, data_pJ, n_trials, selected_individual)

    with st.expander(f'Graphical representation of input'):
        cols = st.columns([K] + [1 for i in range(K)])
        with cols[0]:
            st.bar_chart(data_q, horizontal=True)
            st.bar_chart(data_p, horizontal=True)

    with st.expander(f'Graphical representation of q (individual admixture)'):
        cols = st.columns([1 for i in data_q_out])        
        i = 0
        for d in data_q_out:
            with cols[i]:
                # st.write(d["data"])
                # st.write(d["extension"])
              	st.write(data_q_out)
              	#st.bar_chart(d, horizontal=True)
                i = i+1

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



