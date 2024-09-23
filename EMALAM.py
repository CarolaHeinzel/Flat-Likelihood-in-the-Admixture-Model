import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

import EMALAM_incl_comments_commandline_K2 as em

st.set_page_config(page_title="EMALAM", page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

st.header("EMALAM")

cont0 = cont1 = cont2 = False

if "uploaded_q_file" not in st.session_state:
    st.session_state.uploaded_q_file = ""
if "uploaded_p_file" not in st.session_state:
    st.session_state.uploaded_p_file = ""
if "poss" not in st.session_state:
    st.session_state.poss = ""

col1, col2 = st.columns([1,1])
st.session_state.uploaded_q_file = col1.file_uploader("STRUCTURE output file for allele frequencies (q)")
st.session_state.uploaded_p_file = col2.file_uploader("STRUCTURE output file for individual ancestries (p)")

options = [f"P{i+1}" for i in range(5)]
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

poss = st.selectbox("Function which is optimized", options, index=None, format_func=format_func, key=None, help="See Documentation", on_change=None, args=None, kwargs=None, placeholder="Choose an option", disabled=False, label_visibility="visible")

if poss == "P1":
    st.st.selectbox("Which individual is maximized?")

if st.session_state.uploaded_p_file != "":
    data_p = pd.read_csv(st.session_state.uploaded_p_file, delimiter = " ", header = None)
    M = data_p.shape[0]    
    K = data_p.shape[1]

if st.session_state.uploaded_q_file != "":
    data_q = pd.read_csv(st.session_state.uploaded_q_file, delimiter = " ", header = None)
    N = data_q.shape[0]

submit = st.button("Submit")
if submit:
    print(st.session_state.uploaded_q_file.name)


    # Remove non-polymorphic loci
    data_p = data_p[~(data_p == 1).all(axis=1)]

    pJ = None
    poss = "P2"
    simi = 1
    k_specific = 0
    n_trials = 10 # Number of different initial values for the minimization function

    res_algo_final = em.algo_final(data_q, data_p, K, poss, simi, k_specific, pJ, n_trials)

    st.write(res_algo_final)

#if st.session_state.uploaded_p_file is not None:
#    p_df = pd.read_csv(st.session_state.uploaded_p_file, delimiter = "\t", header = None)
#    M = p_df.shape[0]    
#    K = p_df.shape[1]

#st.session_state.uploaded_q_file = q_file
#if st.session_state.uploaded_q_file is not None:
#    q_df = pd.read_csv(st.session_state.uploaded_q_file, delimiter = "\t", header = None)
#    N = q_df.shape[0]

#cont0 = st.session_state.uploaded_p_file != "" and st.session_state.uploaded_q_file != ""






if cont0:
    if q_df.shape[1] != K:
        st.write("Number of columns in q-file must equal number of columns in p-file (number of populations)")
    else:
        cont1 = True

if cont1:
    cols = st.columns([K] + [1 for i in range(K)])
    with cols[0]:
        st.bar_chart(q_df, horizontal=True)
    for i in range(K):
        with cols[i+1]:
            p_df_loc = pd.DataFrame({'0' : p_df[i], '1': 1-p_df[i]})
            st.bar_chart(p_df_loc, horizontal=True)




# with st.expander(f'Graphical output for structure outputs', expanded = True if st.session_state.uploaded_file is not None else False):


