from pyopenms import AASequence


def parse_sequon(sequon):
    residues_pos, residues_neg = set(), set()
    sequons_pos_neg = sequon.split("!")

    residues_to_exclude = []
    residues_to_include = sequons_pos_neg[0]
    if len(sequons_pos_neg) > 1:
        residues_to_exclude = sequons_pos_neg[1]

    for s in residues_to_include:
        if s == "X":
            residues_pos.update(list('ARNDCEQGHILKMFPSTWYV'))
        elif s in 'ARNDCEQGHILKMFPSTWYV':
            residues_pos.add(s)

    for s in residues_to_exclude:
        if s == "X":
            residues_neg.update(list('ARNDCEQGHILKMFPSTWYV'))
        elif s in 'ARNDCEQGHILKMFPSTWYV':
            residues_neg.add(s)

    return residues_pos - residues_neg


def check_pre_sequon(modified_peptide: AASequence, residue_position, sequons):
    if len(sequons) > residue_position:
        return False

    for i, valid_residues in enumerate(sequons[::-1], 1):
        if modified_peptide.getResidue(residue_position - i).getOneLetterCode() not in valid_residues:
            return False
    return True


def check_post_sequon(modified_peptide: AASequence, residue_position, sequons):
    if residue_position + len(sequons) >= modified_peptide.size():
        return False

    for i, valid_residues in enumerate(sequons, 1):
        if modified_peptide.getResidue(residue_position + i).getOneLetterCode() not in valid_residues:
            return False
    return True


def apply_sequons(modified_peptide: AASequence, sequons):
    if len(sequons) == 0 or modified_peptide.isModified() is False:
        return True
    for i in range(modified_peptide.size()):
        residue = modified_peptide.getResidue(i)
        if residue.isModified():
            mod_name = residue.getModification().getFullId()
            if mod_name in sequons:
                pre_sequon_valid = check_pre_sequon(modified_peptide, i, sequons[mod_name]['pre'])
                post_sequon_valid = check_post_sequon(modified_peptide, i, sequons[mod_name]['post'])

                if pre_sequon_valid is False or post_sequon_valid is False:
                    return False
    return True