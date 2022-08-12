from dataclasses import dataclass
from enum import Enum
from typing import List, Any, Dict
import streamlit as st

from .help_messages import *
from .constants import *
from .sequon_utils import parse_sequon


class TermSpecificity(Enum):
    ANYWHERE = 0
    C_TERM = 1
    N_TERM = 2
    PROTEIN_C_TERM = 3
    PROTEIN_N_TERM = 4


@dataclass
class ModificationParams:
    static_mods: List[str]
    static_peptide_cterm_mods: List[str]
    static_peptide_nterm_mods: List[str]
    static_protein_cterm_mods: List[str]
    static_protein_nterm_mods: List[str]
    variable_mods: List[str]
    variable_peptide_cterm_mods: List[str]
    variable_peptide_nterm_mods: List[str]
    variable_protein_cterm_mods: List[str]
    variable_protein_nterm_mods: List[str]
    max_variable_mods: int
    sequons: Dict[str, Any]

def get_modification_params(mod_db) -> ModificationParams:
    mods = []
    mod_db.getAllSearchModifications(mods)
    mods_by_specificity = {TermSpecificity.ANYWHERE.value: [], TermSpecificity.C_TERM.value: [],
                           TermSpecificity.N_TERM.value: [], TermSpecificity.PROTEIN_C_TERM.value: [],
                           TermSpecificity.PROTEIN_N_TERM.value: []}

    for mod in mods:
        mods_by_specificity[mod_db.getModification(mod).getTermSpecificity()].append(mod)

    add_mass_change_to_mod = lambda l: [f"{s} {mod_db.getModification(s).getDiffMonoMass()}" for s in l]
    convert_mod_strs_into_mods = lambda mod_strs: [mod_db.getModification(" ".join(mod_str.split(" ")[:-1])) for mod_str
                                                   in mod_strs]
    decode_strings = lambda l: [s.decode('utf-8') for s in l]

    with st.expander("Static Modifications"):
        static_mods = st.multiselect('Residue Modifications:',
                                     add_mass_change_to_mod(
                                         decode_strings(mods_by_specificity[TermSpecificity.ANYWHERE.value])),
                                     default=[], help=STATIC_MOD_HELP_COLUMN, key="static_mods")
        static_peptide_cterm_mods = st.multiselect('Peptide C-Term Modifications:',
                                                   add_mass_change_to_mod(
                                                       decode_strings(mods_by_specificity[TermSpecificity.C_TERM.value])),
                                                   default=[], help=STATIC_MOD_HELP_COLUMN, key="static_peptide_cterm_mods")
        static_peptide_nterm_mods = st.multiselect('Peptide N-Term Modifications:',
                                                   add_mass_change_to_mod(
                                                       decode_strings(mods_by_specificity[TermSpecificity.N_TERM.value])),
                                                   default=[], help=STATIC_MOD_HELP_COLUMN, key="static_peptide_nterm_mods")
        static_protein_cterm_mods = st.multiselect('Protein C-Term Modifications:',
                                                   add_mass_change_to_mod(decode_strings(
                                                       mods_by_specificity[TermSpecificity.PROTEIN_C_TERM.value])),
                                                   default=[], help=STATIC_MOD_HELP_COLUMN, key="static_protein_cterm_mods",
                                                   disabled=True)
        static_protein_nterm_mods = st.multiselect('Protein N-Term Modifications:',
                                                   add_mass_change_to_mod(decode_strings(
                                                       mods_by_specificity[TermSpecificity.PROTEIN_N_TERM.value])),
                                                   default=[], help=STATIC_MOD_HELP_COLUMN, key="static_protein_nterm_mods",
                                                   disabled=True)
    with st.expander("Variable Modifications"):
        max_variable_mods = st.number_input('Max variable modifications:', min_value=MIN_VARIABLE_MODIFICATIONS,
                                            max_value=MAX_VARIABLE_MODIFICATIONS, value=DEFAULT_VARIABLE_MODIFICATIONS,
                                            help=MAX_VARIABLE_MODIFICATIONS_HELP_MESSAGE, key="max_var_mods")
        variable_mods = st.multiselect('Residue Modifications:',
                                       add_mass_change_to_mod(
                                           decode_strings(mods_by_specificity[TermSpecificity.ANYWHERE.value])),
                                       default=[], help=VARIABLE_MOD_HELP_COLUMN, key="variable_mods")
        variable_peptide_cterm_mods = st.multiselect('Peptide C-Term Modifications:',
                                                     add_mass_change_to_mod(
                                                         decode_strings(mods_by_specificity[TermSpecificity.C_TERM.value])),
                                                     default=[], help=VARIABLE_MOD_HELP_COLUMN,
                                                     key="variable_peptide_cterm_mods")
        variable_peptide_nterm_mods = st.multiselect('Peptide N-Term Modifications:',
                                                     add_mass_change_to_mod(
                                                         decode_strings(mods_by_specificity[TermSpecificity.N_TERM.value])),
                                                     default=[], help=VARIABLE_MOD_HELP_COLUMN,
                                                     key="variable_peptide_nterm_mods")
        variable_protein_cterm_mods = st.multiselect('Protein C-Term Modifications:',
                                                     add_mass_change_to_mod(decode_strings(
                                                         mods_by_specificity[TermSpecificity.PROTEIN_C_TERM.value])),
                                                     default=[], help=VARIABLE_MOD_HELP_COLUMN,
                                                     key="variable_protein_cterm_mods",
                                                     disabled=True)
        variable_protein_nterm_mods = st.multiselect('Protein N-Term Modifications:',
                                                     add_mass_change_to_mod(decode_strings(
                                                         mods_by_specificity[TermSpecificity.PROTEIN_N_TERM.value])),
                                                     default=[], help=VARIABLE_MOD_HELP_COLUMN,
                                                     key="variable_protein_nterm_mods",
                                                     disabled=True)

        sequons = {}
        if st.checkbox("sequon", help=SEQUON_HELP_MESSAGE):
            mods = convert_mod_strs_into_mods(variable_mods)
            for i, mod in enumerate(mods):
                st.write(mod.getFullId())
                col1, col2 = st.columns(2)
                pre = col1.text_input("Pre residues", key=f"Pre residues{i}")
                post = col2.text_input("Post residues", key=f"Pre residues{i}")
                sequons[mod.getFullId()] = {
                    'pre': [parse_sequon(sequon) for sequon in pre.split("|") if len(parse_sequon(sequon)) > 0],
                    'post': [parse_sequon(sequon) for sequon in post.split("|") if len(parse_sequon(sequon)) > 0]}

    return ModificationParams(static_mods=static_mods,
                              static_peptide_cterm_mods=static_peptide_cterm_mods,
                              static_peptide_nterm_mods=static_peptide_nterm_mods,
                              static_protein_cterm_mods=static_protein_cterm_mods,
                              static_protein_nterm_mods=static_protein_nterm_mods,
                              variable_mods=variable_mods,
                              variable_peptide_cterm_mods=variable_peptide_cterm_mods,
                              variable_peptide_nterm_mods=variable_peptide_nterm_mods,
                              variable_protein_cterm_mods=variable_protein_cterm_mods,
                              variable_protein_nterm_mods=variable_protein_nterm_mods,
                              max_variable_mods=max_variable_mods,
                              sequons=sequons)