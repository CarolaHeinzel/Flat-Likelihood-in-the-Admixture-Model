import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

import EMALAM_incl_comments_commandline_K2 as em

st.header("EMALAM")

col1, col2 = st.columns([1,1]) 
p_file = col1.file_uploader("STRUCTURE output file for allele frequencies (p)")
q_file = col2.file_uploader("STRUCTURE output file for individual ancestries (q)")

st.session_state.uploaded_p_file = p_file
if st.session_state.uploaded_p_file is not None:
    p_df = pd.read_csv(st.session_state.uploaded_p_file, delimiter = " ", header = None)
    p_fig = px.pie(p_df)
    fig = st.plotly_chart(p_fig, theme=None)



st.session_state.uploaded_q_file = q_file
if st.session_state.uploaded_q_file is not None:
    q_df = pd.read_csv(st.session_state.uploaded_q_file, delimiter = " ", header = None)
    # col2.write(q_df)
    # st.write(tuple(q_df[0]))
    st.bar_chart(q_df, horizontal=True)

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

submit = st.button("Submit")
# if submit:


# with st.expander(f'Graphical output for structure outputs', expanded = True if st.session_state.uploaded_file is not None else False):


