from dataclasses import dataclass

import streamlit as st
from pyopenms import ProteaseDB

from .constants import *
from.help_messages import *


@dataclass
class DigestionParams:
    missed_cleavages: int
    min_peptide_length: int
    max_peptide_length: int
    protease: str
    specificity: str
    decoy_flag: str


def get_supported_proteases():
    names = []
    ProteaseDB().getAllNames(names)
    return [str(name.decode("utf-8")) for name in names]


def get_digestion_params() -> DigestionParams:
    missed_cleavages = int(st.number_input('Missed Cleavages:', min_value=MIN_MISSED_CLEAVAGES,
                                           max_value=MAX_MISSED_CLEAVAGES, value=DEFAULT_MISSED_CLEAVAGES,
                                           help=MISSED_CLEAVAGES_HELP_MESSAGE))
    min_peptide_length = int(st.number_input('Min Peptide Length:', min_value=MIN_PEPTIDE_LENGTH,
                                             max_value=MAX_PEPTIDE_LENGTH, value=DEFAULT_MIN_PEPTIDE_LENGTH,
                                             help=MIN_PEPTIDE_LENGTH_HELP_MESSAGE))
    max_peptide_length = int(st.number_input('Max Peptide Length:', min_value=MIN_PEPTIDE_LENGTH,
                                             max_value=MAX_PEPTIDE_LENGTH, value=DEFAULT_MAX_PEPTIDE_LENGTH,
                                             help=MAX_PEPTIDE_LENGTH_HELP_MESSAGE))
    supported_proteases = get_supported_proteases()
    protease = st.selectbox('Protease:', supported_proteases, index=supported_proteases.index('Trypsin'),
                            help=PROTEASE_HELP_MESSAGE)
    specificity = st.selectbox('Protease Specificity:', ('None', 'semi', 'full'), index=2,
                               help=SPECIFICITY_HELP_MESSAGE)

    return DigestionParams(missed_cleavages=missed_cleavages,
                          min_peptide_length=min_peptide_length,
                          max_peptide_length=max_peptide_length,
                          protease=protease,
                          specificity=specificity,
                          decoy_flag="DECOY_")