from pyopenms import AASequence


def get_mod_string(sequence):
    p = AASequence.fromString(sequence)
    mods = []
    if p.hasNTerminalModification():
        mods.append("0")
        mods.append(p.getNTerminalModificationName())

    for i in range(p.size()):
        residue = p.getResidue(i)
        mod = residue.getModificationName()
        if mod:
            mods.append(str(i))
            mods.append(mod)

    if p.hasCTerminalModification():
        mods.append("-1")
        mods.append(p.getCTerminalModificationName())

    return "|".join(mods)


def get_mod_and_locations(sequence):
    p = AASequence.fromString(sequence)
    mods = []
    locations = []
    if p.hasNTerminalModification():
        locations.append("0")
        mods.append(f"{p.getNTerminalModificationName()}@{p.getResidue(0).getOneLetterCode()}")

    for i in range(p.size()):
        residue = p.getResidue(i)
        mod = residue.getModificationName()
        if mod:
            locations.append(str(i))
            mods.append(f"{mod}@{residue.getOneLetterCode()}")

    if p.hasCTerminalModification():
        locations.append("-1")
        mods.append(f"{p.getCTerminalModificationName()}@{p.getResidue(p.size()-1).getOneLetterCode()()}")

    return mods, locations