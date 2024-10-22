import streamlit as st
import numpy as np
import altair as alt
import pandas as pd
import emalam as em

def plot_q(q_df):
    st.altair_chart(alt.Chart(q_df).mark_bar().encode(
        x=alt.X("ind:N", sort=None, title=""),
        y=alt.Y("ia:Q", sort = 'descending', title="", scale=alt.Scale(domain=[0, 1])),
        color=alt.Color('pop:N', scale=alt.Scale(scheme='category10')),
        tooltip=['pop:N', 'ia:Q']  
    ).properties(), use_container_width=True)                    

def plot_p(p_df):
    K = em.get_K_from_hatp_dict(hatp_dict)
    for k in range(K):
        st.write(f"Population {k}")
        st.altair_chart(alt.Chart(p_df[p_df["pop"] == k]).mark_bar().encode(
            x=alt.X("id:N", title=""),
            y=alt.Y("p:Q", sort = 'descending', title=""),
            color=alt.Color('allele:N', scale=alt.Scale(scheme='category10')),
            tooltip=['marker:N', 'allele:N', 'p:Q']  
        ).properties(), use_container_width=True)

default_structure_path_K3 = 'Example_Input/CEU_IBS_TSI_enhanced_corr_K3_f'
default_structure_path_K3_url = "https://github.com/CarolaHeinzel/Flat-Likelihood-in-the-Admixture-Model/blob/main/Example_Input/CEU_IBS_TSI_enhanced_corr_K3_f"
default_structure_path_K2 = 'Example_Input/CEU_IBS_TSI_enhanced_corr_f'
default_structure_path_K2_url = "https://github.com/CarolaHeinzel/Flat-Likelihood-in-the-Admixture-Model/blob/main/Example_Input/CEU_IBS_TSI_enhanced_corr_f"

default_q_path = 'Example_Input/q_CEU_IBS_TSI_K3'
default_p_path = 'Example_Input/p_CEU_IBS_TSI_K3'
default_q_path_url = "https://github.com/CarolaHeinzel/Flat-Likelihood-in-the-Admixture-Model/blob/main/Example_Input/q_CEU_IBS_TSI_K3"
default_p_path_url = "https://github.com/CarolaHeinzel/Flat-Likelihood-in-the-Admixture-Model/blob/main/Example_Input/p_CEU_IBS_TSI_K3"

n = 10 # number of trials for optimization

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
    st.session_state.inds = []
if "ids" not in st.session_state:
    st.session_state.ids = None

lines = None

st.header("EMALAM")
st.write("See [biorxiv](https://www.biorxiv.org/content/10.1101/2024.10.18.619150v1) for a reference what this webpage is about.")

################################
## Choose which input is used ##
################################

# At the end of this section, we have hatq and hatp

default_options = [
    f"Example with K=3 [STRUCTURE output data]({default_structure_path_K3_url})",
    f"Example with K=2 [STRUCTURE output data]({default_structure_path_K2_url})",
    f"Example ADMIXTURE output data [Q]({default_q_path_url}) and [P]({default_p_path_url})", 
    "Upload data"]
default = st.radio(
    "Which file(s) should be used?",
    options = default_options,
    index = 0,
)

if default == default_options[0]:
    use_structure_file = True
    lines = em.load_structure_file(default_structure_path_K3)

if default == default_options[1]:
    use_structure_file = True
    lines = em.load_structure_file(default_structure_path_K2)
if default == default_options[2]:
    use_structure_file = False
    st.session_state.uploaded_q_file = default_q_path
    st.session_state.uploaded_p_file = default_p_path
if default == default_options[3]:
    data_type_options = ["STRUCTURE output (single file)", "Admixture output (Q and P files)"]
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
        col0, col1 = st.columns([1,1])
        st.session_state.uploaded_q_file = col0.file_uploader("ADMIXTURE output file for individual ancestries (Q)")
        st.session_state.uploaded_p_file = col1.file_uploader("ADMIXTURE output file for allele frequencies (P)")

cont = (use_structure_file and lines is not None) or (not use_structure_file and st.session_state.uploaded_q_file is not None and st.session_state.uploaded_p_file is not None)

if cont:
    if use_structure_file:
        hatq_dict = em.get_q(lines)
        hatq_df = em.to_df(hatq_dict)
        st.session_state.hatq = np.array(hatq_df)
        st.session_state.ind_ids = hatq_dict.keys()
        hatp_dict = em.get_p(lines)    
    else:
        st.session_state.hatq = np.array(em.load_q_file(st.session_state.uploaded_q_file))
        st.session_state.ind_ids = range(st.session_state.hatq.shape[0])
        hatp_dict = em.load_p_file(st.session_state.uploaded_p_file)

    hatp_df = em.to_df(hatp_dict)
    hatp = np.array(hatp_df)
    st.session_state.hatp = np.array(hatp_df)

############################
## Choose target function ##
############################

if cont: #lines is not None:
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

    if target == target_options[1]:
        pop_options = range(st.session_state.hatq.shape[1])
        st.session_state.popl = st.selectbox("...for which population i", 
                        options = pop_options,
                        index = 0,
                        help = "i is the column number in the structure output file, starting with 0")

    showq = st.toggle("Show me the individual admixtures for all individuals", True)
    showp = st.toggle("Show me the allele frequencies at all markers", False)
    notall = st.toggle("Minimize/maximize target function only for a subset of individuals", False)
    all = not notall
    
    if all:
        st.session_state.inds = [True for id in st.session_state.ind_ids]
    else:
        st.session_state.ids = st.multiselect(
            "Average the target function on the following individuals",
            options = st.session_state.ind_ids, placeholder = "Please choose one or more individuals.")
        st.session_state.inds = [id in st.session_state.ids for id in st.session_state.ind_ids]
            
###################
## Run optimizer ##    
###################

if cont: # lines is not None:
    hatq = st.session_state.hatq
    hatp = st.session_state.hatp
    inds = st.session_state.inds
    popl = st.session_state.popl
    ind_ids = st.session_state.ind_ids
    ids = st.session_state.ids
    K = hatq.shape[1]
    N = hatq.shape[0]
    M = hatp.shape[0]

    submit = st.button("Submit")
    if submit:
        if st.session_state.ids == [] and notall == True:
             notall = False
        with st.spinner('The computer is calculating...'):
            if target == target_options[0]:
                S_min = em.find_S(em.mean_entropy, hatq, hatp, n, (hatq, inds), em.mean_entropy_jac)
                S_max = em.find_S(em.neg_mean_entropy, hatq, hatp, n, (hatq, inds), em.neg_mean_entropy_jac) 
                #q_min = em.find_q(em.mean_entropy, hatq, hatp, n, (hatq, inds), em.mean_entropy_jac)
                #q_max = em.find_q(em.neg_mean_entropy, hatq, hatp, n, (hatq, inds), em.neg_mean_entropy_jac) 
            else:
                S_min = em.find_S(em.mean_size, hatq, hatp, n, (hatq, popl, inds), em.mean_size_jac)
                S_max = em.find_S(em.neg_mean_size, hatq, hatp, n, (hatq, popl, inds), em.neg_mean_size_jac)
                #q_min = em.find_q(em.mean_size, hatq, hatp, n, (hatq, popl, inds), em.mean_size_jac)
                #q_max = em.find_q(em.neg_mean_size, hatq, hatp, n, (hatq, popl, inds), em.neg_mean_size_jac)
            q_min = hatq.dot(S_min)
            q_max = hatq.dot(S_max)

        ####################
        ## Output results ##    
        ####################

        with st.expander("### Range of individual ancestries"):
            if notall:
                st.write("Individual admixtures of selected individuals")
                for i in range(len(ids)):
                    col0, col1 = st.columns([1,4])
                    with col0:
                        st.write(f"Individual {ids[i]}")
                    with col1:
                        q = np.array([hatq[i], q_min[i], q_max[i]]).flatten()
                        q_df = pd.DataFrame({
                                "ia": q,
                                "cat": ["start" for i in range(K)] + ["min" for i in range(K)] + ["max" for i in range(K)], 
                                "pop": [ j for i in range(3) for j in range(K)]
                            })
        #                    st.write(q_df)
                        st.altair_chart(alt.Chart(q_df).mark_bar().encode(
                            x=alt.Y("ia:Q", title="", scale=alt.Scale(domain=[0, 1])),
                            y=alt.X("cat:N", sort=None, title=""),
                            color=alt.Color('pop:N', scale=alt.Scale(scheme='category10')),
                            tooltip=['pop:N', 'ia:Q']  
                        ).properties(width=700), use_container_width=False)


            else:
                st.write("## Individual admixtures of all individuals")
                st.write(f"N = {N}, K = {K}, M = {M}")
                st.write("#### Initial values")
                col0, col1 = st.columns([10,2])
                q_df, hatq_pivot = em.get_q_for_plot(hatq)
                with col0:
                    plot_q(q_df)
                with col1:
                    st.write(hatq_pivot)
                st.write("#### Minimum")
                col0, col1 = st.columns([10,2])
                q_df, hatq_pivot = em.get_q_for_plot(q_min)
                with col0:
                    plot_q(q_df)
                with col1:
                    st.write(hatq_pivot)
                st.write("#### Maximum")
                col0, col1 = st.columns([10,2])
                q_df, hatq_pivot = em.get_q_for_plot(q_max)
                with col0:
                    plot_q(q_df)
                with col1:
                    st.write(hatq_pivot)

        if showp:
            with st.expander("Allele frequencies at all markers"):
                T_min = np.linalg.pinv(S_min)
                T_max = np.linalg.pinv(S_max)
                p_min_dict = { key: value.dot(T_min.T) for key, value in hatp_dict.items() } 
                p_max_dict = { key: value.dot(T_max.T) for key, value in hatp_dict.items() } 

                st.write("#### Initial values")
                col0, col1 = st.columns([10,2])
                p_df, hatp_pivot = em.get_p_for_plot(hatp_dict)
                with col0:
                    plot_p(p_df)
                with col1:
                    st.write(hatp_pivot)
                st.write("#### p for minimal target function")
                col0, col1 = st.columns([10,2])
                p_df, hatp_pivot = em.get_p_for_plot(p_min_dict)
                with col0:
                    plot_p(p_df)
                with col1:
                    st.write(hatp_pivot)
                st.write("#### p for maximal target function")
                p_df, hatp_pivot = em.get_p_for_plot(p_max_dict)
                col0, col1 = st.columns([10,2])
                with col0:
                    plot_p(p_df)
                with col1:
                    st.write(hatp_pivot)


