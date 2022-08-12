from functools import partial

import pandas as pd
import streamlit as st

from pyopenms import *

from modification.utils import get_modified_peptides
from modification.params import get_modification_params
from modification.sequon_utils import apply_sequons

mod_db = ModificationsDB()

st.header("Visualize Peptide Modifications")
st.write("Can help to understand the complexity which modifications add to database searches.")

st.subheader("Input")
peptide = st.text_input('Peptide Sequence', 'ARNDCEQGHILKMFPSTWYV')

st.subheader("Params")
modification_params = get_modification_params(mod_db)

convert_mod_strs_into_mods = lambda mod_strs: [mod_db.getModification(" ".join(mod_str.split(" ")[:-1])) for mod_str
                                               in mod_strs]
get_modified_peptides_partial = partial(get_modified_peptides,
                                        static_mods=convert_mod_strs_into_mods(modification_params.static_mods),
                                        static_peptide_cterm_mods=
                                        convert_mod_strs_into_mods(modification_params.static_peptide_cterm_mods),
                                        static_peptide_nterm_mods=
                                        convert_mod_strs_into_mods(modification_params.static_peptide_nterm_mods),
                                        variable_mods=convert_mod_strs_into_mods(modification_params.variable_mods),
                                        variable_peptide_cterm_mods=
                                        convert_mod_strs_into_mods(modification_params.variable_peptide_cterm_mods),
                                        variable_peptide_nterm_mods=
                                        convert_mod_strs_into_mods(modification_params.variable_peptide_nterm_mods),
                                        max_variable_mods=modification_params.max_variable_mods)
st.subheader("Output")
peptide = AASequence.fromString(peptide)

modified_peptides = get_modified_peptides_partial(peptide)
modified_peptides = [p for p in modified_peptides if apply_sequons(p, modification_params.sequons) is True]
data = {'peptide': [str(p) for p in modified_peptides]}
mod_df = pd.DataFrame(data)
st.dataframe(mod_df)
