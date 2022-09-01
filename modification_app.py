import pandas as pd
import streamlit as st

from pyopenms import *

from modification.params import get_modification_params
from modification.sequon_utils import apply_sequons
from modification.utils import get_modified_peptides_partial

mod_db = ModificationsDB()

st.header("Visualize Peptide Modifications")
with st.expander("Help"):
    st.write("Can help to understand the complexity which modifications add to database searches.")

st.subheader("Input")
peptide = st.text_input('Peptide Sequence', 'ARNDCEQGHILKMFPSTWYV')

st.subheader("Params")
modification_params = get_modification_params(mod_db)

get_modified_peptides = get_modified_peptides_partial(mod_db, modification_params)
st.subheader("Output")
peptide = AASequence.fromString(peptide)

modified_peptides = get_modified_peptides(peptide)
modified_peptides = [p for p in modified_peptides if apply_sequons(p, modification_params.sequons) is True]
data = {'peptide': [str(p) for p in modified_peptides]}
mod_df = pd.DataFrame(data)
st.table(mod_df)
