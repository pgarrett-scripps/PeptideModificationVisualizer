from functools import partial

from pyopenms import AASequence


def apply_residue_static_mods(static_mods, peptide):
    for mod in static_mods:
        for i in range(peptide.size()):
            residue = peptide.getResidue(i)
            if residue.getOneLetterCode() == mod.getOrigin():
                peptide.setModification(i, mod.getFullId())


def apply_static_c_term_mods(mods, peptide):
    term_residue = peptide.getResidue(peptide.size() - 1).getOneLetterCode()
    for mod in mods:
        origin_residue = mod.getOrigin()
        if origin_residue == 'X' or origin_residue == term_residue:
            peptide.setCTerminalModification(mod.getFullId())


def apply_static_n_term_mods(mods, peptide):
    term_residue = peptide.getResidue(0).getOneLetterCode()
    for mod in mods:
        origin_residue = mod.getOrigin()
        if origin_residue == 'X' or origin_residue == term_residue:
            peptide.setNTerminalModification(mod.getFullId())


def apply_static_mods(peptide, residue_mods, cterm_mods, nterm_mods):
    apply_residue_static_mods(residue_mods, peptide)
    apply_static_c_term_mods(cterm_mods, peptide)
    apply_static_n_term_mods(nterm_mods, peptide)


def rec_mod_builder(mods, peptide, i, mod_count, max_mod_count, results):
    if len(mods) == 0:  # skip for speed up
        return results[peptide]

    # print(str(peptide), i, mod_count, max_mod_count, results)
    if i == peptide.size() - 1:
        results.append(peptide)
        return

    residue = peptide.getResidue(i)
    for mod in mods:

        if mod_count >= max_mod_count:
            break

        if residue.getOneLetterCode() == mod.getOrigin():
            new_peptide = AASequence(peptide)
            new_peptide.setModification(i, mod.getFullName())
            rec_mod_builder(mods, new_peptide, i + 1, mod_count + 1, max_mod_count, results)

    rec_mod_builder(mods, peptide, i + 1, mod_count, max_mod_count, results)


def apply_variable_c_term_mods(mods, peptide, results):
    term_residue = peptide.getResidue(peptide.size() - 1).getOneLetterCode()
    for mod in mods:
        origin_residue = mod.getOrigin()
        if origin_residue == 'X' or origin_residue == term_residue:
            new_peptide = AASequence(peptide)
            new_peptide.setCTerminalModification(mod.getFullId())
            results.append(new_peptide)


def apply_variable_n_term_mods(mods, peptide, results):
    term_residue = peptide.getResidue(0).getOneLetterCode()
    for mod in mods:
        origin_residue = mod.getOrigin()
        if origin_residue == 'X' or origin_residue == term_residue:
            new_peptide = AASequence(peptide)
            new_peptide.setNTerminalModification(mod.getFullId())
            results.append(new_peptide)


def apply_variable_mods(peptide, residue_mods, cterm_mods, nterm_mods, max_mod_count):
    c_term_var_peptides = [peptide]
    apply_variable_c_term_mods(cterm_mods, peptide, c_term_var_peptides)

    n_term_var_peptides = [peptide]
    apply_variable_n_term_mods(nterm_mods, peptide, n_term_var_peptides)

    term_peptides = []
    for c_term_var_peptide in c_term_var_peptides:
        for n_term_var_peptide in n_term_var_peptides:
            combined_peptide = AASequence(c_term_var_peptide)
            combined_peptide.setNTerminalModification(n_term_var_peptide.getNTerminalModificationName())
            term_peptides.append(combined_peptide)

    results = []
    for p in term_peptides:
        rec_mod_builder(residue_mods, p, 0, 0, max_mod_count, results)

    return results


def get_modified_peptides(peptide, static_mods, static_peptide_cterm_mods, static_peptide_nterm_mods, variable_mods,
                          variable_peptide_cterm_mods, variable_peptide_nterm_mods, max_variable_mods):
    if len(static_mods) != 0 or len(static_peptide_cterm_mods) != 0 or len(static_peptide_nterm_mods) != 0:
        apply_static_mods(peptide, static_mods, static_peptide_cterm_mods, static_peptide_nterm_mods)

    if len(variable_mods) == 0 and len(variable_peptide_cterm_mods) == 0 and len(variable_peptide_nterm_mods) == 0:
        return [peptide]

    modified_peptides = apply_variable_mods(peptide, variable_mods, variable_peptide_cterm_mods,
                                            variable_peptide_nterm_mods, max_variable_mods)
    return modified_peptides

def get_modified_peptides_partial(mod_db, modification_params):
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
    return get_modified_peptides_partial
