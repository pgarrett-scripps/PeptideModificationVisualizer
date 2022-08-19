from pathlib import Path

import streamlit as st
from senpy.dtaSelectFilter.parser2 import read_file


def map_locus_to_sequnce_from_fasta(fasta_lines):
    locus_to_sequence_map = {}
    locus = None
    for line in fasta_lines:
        if line == "":
            continue
        elif line[0] == ">":  # new protein
            locus = line.rstrip().split(" ")[0].replace(">", "")
            description = " ".join(line.rstrip().split(" ")[1:])
            locus_to_sequence_map[locus] = {'sequence': "", 'description': description}
        else:  # protein sequence
            locus_to_sequence_map[locus]['sequence'] += line.rstrip()
    return locus_to_sequence_map


def fasta_from_locus_to_sequnce_map(locus_to_sequnce_map):
    lines = []
    for locus in locus_to_sequnce_map:
        lines.append(f'>{locus} {locus_to_sequnce_map[locus]["description"]}\n')
        lines.append(f'{locus_to_sequnce_map[locus]["sequence"]}\n')
    return lines


st.title("FASTA Filter")
st.write("Generates a new Fasta file from filtered results")

fasta_file = st.file_uploader("FASTA", ".fasta")
dta_filter_files = st.file_uploader("DTASelect-filter.txt", ".txt", accept_multiple_files=True)

if fasta_file and dta_filter_files:

    fasta_lines = fasta_file.getvalue().decode("utf-8").split("\n")
    results_list = []
    for dta_filter_file in dta_filter_files:
        filter_lines = dta_filter_file.getvalue().decode("utf-8").split("\n")
        _, results, _ = read_file(filter_lines)
        results_list.extend(results)

    locus_to_sequence_map = map_locus_to_sequnce_from_fasta(fasta_lines)
    protein_locuses = list({result.protein.locus for result in results_list})
    dta_filter_locus_to_sequence_map = {locus: locus_to_sequence_map[locus] for locus in protein_locuses}

    new_fasta_lines = fasta_from_locus_to_sequnce_map(dta_filter_locus_to_sequence_map)
    st.download_button("Download FASTA", "".join(new_fasta_lines), file_name=f"{Path(fasta_file.name).stem}_filter.fasta")
