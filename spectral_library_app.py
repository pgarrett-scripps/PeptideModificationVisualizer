import time
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from pyopenms import ModificationsDB, ProteaseDigestion, AASequence

from digestion.params import get_digestion_params
from fasta.util import load_fasta
from modification.utils import get_modified_peptides_partial
from modification.params import get_modification_params
from modification.sequon_utils import apply_sequons
from machine_learning.utils import get_mod_and_locations
from peptdeep.pretrained_models import ModelManager

def flatten(l):
    return [item for sublist in l for item in sublist]

def update_progress_bar(progress, last_progress, total_items, update_frequency=0.1):
    return progress - last_progress >= total_items * update_frequency

def convert_df(df_to_download):
    return df_to_download.to_csv(index=False).encode('ascii')

mod_db = ModificationsDB()

fasta_file = st.file_uploader("Choose a fasta file", type=".fasta")


with st.expander("Digestion"):
    digestion_params = get_digestion_params()

modification_params = get_modification_params(mod_db)

min_charge, max_charge = 1,2
instrument_options = ["Lumos","QE","timsTOF","SciexTOF"]
with st.expander("Deep Learning"):
    retention_time_prediction = st.checkbox('Retention Time', value=False)
    ccs_prediction = st.checkbox('CCS/OOK0', value=False)
    fragment_ion_prediction = st.checkbox('Fragment Ion Intensities', value=False)
    if fragment_ion_prediction:
        instrument = st.selectbox("Instrument", options=instrument_options, index=0)
        nce = st.number_input("NCE", min_value=0.0, max_value=1.0, value=0.3)

    if ccs_prediction or fragment_ion_prediction:
        min_charge, max_charge = st.slider("min/max charge", min_value=1, max_value=10, value=(2, 4))


if st.button("Generate Spectral Library"):
    if not fasta_file:
        st.warning('Upload a FASTA file!')
    else:
        st.warning("Don't change params or close tab else the current process will be terminated.")

        dig = ProteaseDigestion()
        dig.setEnzyme(digestion_params.protease)
        dig.setMissedCleavages(digestion_params.missed_cleavages)
        dig.setSpecificity({'None': 0, 'semi': 1, 'full': 2}[digestion_params.specificity])
        get_modified_peptides = get_modified_peptides_partial(mod_db, modification_params)

        state = st.text("Staring Spectral Library Generation...")
        my_bar = st.progress(0)

        st.subheader("Fasta Metrics:")
        protein_col, peptide_col, unique_peptide_col = st.columns(3)

        state.text('Reading FASTA...')
        start_time = time.time()
        fasta_entries = load_fasta(fasta_file)

        protein_col.metric(label="proteins", value=len(fasta_entries))

        state.text("Digesting Proteins...")
        my_bar.progress(0)
        last_progress_update, total_proteins = 0, len(fasta_entries)
        peptide_count = 0
        sequence_peptide_map, sequence_to_protein_map = {}, {}
        for i, protein in enumerate(fasta_entries):
            bsa = AASequence.fromString(protein.sequence)
            peptides = []
            dig.digest(bsa, peptides, digestion_params.min_peptide_length, digestion_params.max_peptide_length)
            peptides = [peptide for peptide in peptides if 'X' not in peptide.toString()]
            modified_peptides = flatten([get_modified_peptides_partial(peptide) for peptide in peptides])
            modified_peptides = [p for p in modified_peptides if apply_sequons(p, modification_params.sequons) is True]
            peptide_count += len(modified_peptides)
            for peptide in modified_peptides:
                peptide_sequence = peptide.toString()
                sequence_peptide_map[peptide_sequence] = peptide
                if peptide_sequence not in sequence_to_protein_map:
                    sequence_to_protein_map[peptide_sequence] = []
                sequence_to_protein_map[peptide_sequence].append(protein)

            if my_bar and update_progress_bar(i, last_progress_update, total_proteins, update_frequency=0.01):
                my_bar.progress(i / total_proteins)
                last_progress_update = i

        peptide_col.metric(label="peptides", value=peptide_count)
        unique_peptide_col.metric(label="unique peptides", value=len(sequence_peptide_map))

        data = {'peptide':[],'protein':[], 'charge':[]}
        for peptide_sequence in sequence_to_protein_map:
            protein_ids = [protein.identifier for protein in sequence_to_protein_map[peptide_sequence]]
            for charge in range(min_charge, max_charge+1):
                data['peptide'].append(peptide_sequence)
                data['protein'].append(" ".join(protein_ids))
                data['charge'].append(charge)

        df = pd.DataFrame(data)
        df['mod_sequence'] = [AASequence.fromString(seq).toUniModString() for seq in df['peptide'].values]
        df['sequence'] = [AASequence.fromString(seq).toUnmodifiedString() for seq in df['mod_sequence'].values]
        mod_locs = [get_mod_and_locations(seq) for seq in df['mod_sequence'].values]
        df['mods'] = [";".join(mod_loc[0]) for mod_loc in mod_locs]
        df['mod_sites'] = [";".join(mod_loc[1]) for mod_loc in mod_locs]

        with st.spinner("Running Machine Learning..."):
            model_mgr = ModelManager()
            model_mgr.load_installed_models()
            if retention_time_prediction is True:
                state.text("Predicting Retention Time...")
                model_mgr.predict_rt(df)
            if ccs_prediction is True:
                state.text("Predicting Mobility...")
                model_mgr.predict_mobility(df)
            if fragment_ion_prediction is True:
                state.text("Predicting Fragment Ion Intensities...")
                df['nce'] = nce
                df['instrument'] = instrument
                df['intensities'] = None
                df_intensities = model_mgr.predict_ms2(df)
                for index, row in df.iterrows():
                    flatten_intensities = np.concatenate(
                        df_intensities.iloc[row.frag_start_idx:row.frag_end_idx].values)
                    df.loc[index, 'intensities'] = np.array2string(flatten_intensities, precision=2, separator=',',
                                                                   suppress_small=True)

        st.dataframe(df)

        st.download_button(label="Download",

                           data=convert_df(df),
                           file_name=f'{Path(fasta_file.name).stem}_spectral_lib.csv',
                           mime='text/csv', )

        st.balloons()
