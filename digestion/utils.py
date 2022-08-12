def get_target_decoy_proteins_target_decoy_fasta(fasta_entries, decoy_flag):
    target_protein_entries, decoy_protein_entries = [], []
    for i, fasta_entry in enumerate(fasta_entries):
        if decoy_flag in fasta_entry.identifier:
            decoy_protein_entries.append(fasta_entry)
        else:
            target_protein_entries.append(fasta_entry)

    return target_protein_entries, decoy_protein_entries