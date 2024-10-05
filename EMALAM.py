import streamlit as st
import plotly.express as px
import matplotlib.pyplot as plt

import tools 

default_structure_path = 'Example_Input/output_structure_K3_f'
n = 10 # number of trials for optimization

st.set_page_config(page_title="EMALAM", page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

# Initialize session state variables
if "uploaded_structure_file" not in st.session_state:
    st.session_state.uploaded_structure_file = None
if "uploaded_q_file" not in st.session_state:
    st.session_state.uploaded_q_file = None
if "uploaded_p_file" not in st.session_state:
    st.session_state.uploaded_p_file = None
if "uploaded_pJ_file" not in st.session_state:
    st.session_state.uploaded_pJ_file = None
if "poss" not in st.session_state:
    st.session_state.poss = "P4" # Initial Value
if "data_uploaded" not in st.session_state:
    st.session_state.data_uploaded = None
if "hatp" not in st.session_state:
    st.session_state.hatp = None
if "hatq" not in st.session_state:
    st.session_state.hatq = None
if "ind_ids" not in st.session_state:
    st.session_state.ind_ids = None
if "popl" not in st.session_state:
    st.session_state.popl = None
if "inds" not in st.session_state:
    st.session_state.inds = None
if "ids" not in st.session_state:
    st.session_state.ids = None

lines = None

default_p_path = 'Example_Input/p_all_K2'
default_q_path = 'Example_Input/q_K2'
default_pJ_path = 'Example_Input/pJ_K2'
default_STRUCTURE_path = 'Example_Input/output_structure_f'

st.header("EMALAM")
st.write("See xxx biorxiv for a reference what this webpage is about.")

################################
## Choose which input is used ##
################################

# At the end of this section, we have hatq and hatp

default_options = ["Example [STRUCTURE output data](github.com)", "Example [p- and q-matrizes](github.com)", "Upload data"]
default = st.radio(
    "Which file(s) should be used?",
    options = default_options,
    index = 0,
)

if default == default_options[0]:
    use_structure_file = True
    lines = tools.load_structure_file(default_structure_path)
if default == default_options[1]:
    use_structure_file = False
    st.write("Sorry, no implementation yet!")
if default == default_options[2]:
    data_type_options = ["STRUCTURE Output", "Other"]
    data_type = st.selectbox(
        "Select the type of input data:",
        options=data_type_options,  # Add other options as needed
        index=0  # Set the default selected option to "STRUCTURE"
    )
    
    # File uploaders based on input type
    if data_type == data_type_options[0]:
        use_structure_file = True
        st.session_state.uploaded_structure_file = st.file_uploader("STRUCTURE output file")
        if st.session_state.uploaded_structure_file:
            lines = st.session_state.uploaded_structure_file.readlines()
        
    else:
        use_structure_file = False
        col1, col2, col3 = st.columns([1,1,1])
        st.session_state.uploaded_q_file = col1.file_uploader("STRUCTURE output file for individual ancestries (q)")
        st.session_state.uploaded_p_file = col2.file_uploader("STRUCTURE output file for allele frequencies (p)")
        st.session_state.uploaded_pJ_file = col3.file_uploader("STRUCTURE output file for allele frequencies (p) for J >= 3")
        st.write("Sorry, no implementation yet!")
        lines = None
    
if lines is not None:
    hatq_df = tools.get_q(lines)
    st.session_state.ind_ids = hatq_df.keys()
    st.session_state.hatq = tools.to_df(hatq_df)    
    hatp_df = tools.get_p(lines)
    st.session_state.hatp = tools.to_df(hatp_df)


############################
## Choose target function ##
############################

if lines is not None:
    st.write("### Selection of target function for optimization")
    st.write("What should be done?")
    
    target_options = [
        "Minimize/Maximize entropy", 
        "Minimize/Maximize contribution from population i"]

    target = st.selectbox(
        "Select target function", 
        options=target_options,  # Add other options as needed
        index=0,  # Set the default selected option to "minimize entropy"
        help = "Maximal entropy favours admixed individuals, minimal entropy favours non-admixed individuals."
    )
    col1, col2 = st.columns([1,1])
    with col1:
        st.session_state.ids = st.multiselect(
            "Average the target function on the following individuals",
            options = st.session_state.ind_ids, placeholder = "Please choose one or more individuals.")
        st.session_state.inds = [id in st.session_state.ids for id in st.session_state.ind_ids]
        
    with col2:
        if target == target_options[1]:
            pop_options = range(st.session_state.hatq.shape[1])
            st.session_state.popl = st.selectbox("...for which population i", 
                               options = pop_options,
                               index = 0,
                               help = "i is the column number in the structure output file, starting with 0")


###################
## Run optimizer ##    
###################

if lines is not None:
    if st.session_state.hatq is not None and st.session_state.hatq is not None and any(st.session_state.inds) == True and (target == target_options[0] or st.session_state.popl is not None):
        submit = st.button("Submit")
        if submit:
            with st.spinner('The computer is calculating...'):
                if target == target_options[0]:
                    fun = (lambda x : tools.mean_entropy(x, st.session_state.hatq, st.session_state.inds))
                    q_min = tools.find_q(fun, st.session_state.hatq, st.session_state.hatp, st.session_state.inds, n)
                    fun = (lambda x : -tools.mean_entropy(x, st.session_state.hatq, st.session_state.inds))
                    q_max = tools.find_q(fun, st.session_state.hatq, st.session_state.hatp, st.session_state.inds, n)
                else:
                    fun = (lambda x : tools.mean_size(x, st.session_state.hatq, st.session_state.popl))
                    q_min = tools.find_q(fun, st.session_state.hatq, st.session_state.hatp, st.session_state.inds, n)
                    fun = (lambda x : -tools.mean_size(x, st.session_state.hatq, st.session_state.popl))
                    q_max = tools.find_q(fun, st.session_state.hatq, st.session_state.hatp, st.session_state.inds, n)
            col0, col1, col2, col3 = st.columns([1,1,1,1])
            with col1:
                st.write("IA (upload)", help = "IA = Individual Ancestry")
            with col2:
                st.write("IA minimized", help = "IA = Individual Ancestry")
            with col3:
                st.write("IA maximized", help = "IA = Individual Ancestry")
            for i in range(len(st.session_state.ids)):
                col0, col1, col2, col3 = st.columns([1,1,1,1])
                with col0:
                    st.write(f"Individual {st.session_state.ids[i]}")
                with col1:
                    hatq = tools.subset_inds(st.session_state.hatq, st.session_state.inds)
                    st.write(hatq[i])
                with col2:
                    st.write(q_min[i])
                with col3:
                    st.write(q_max[i])


