SEQUON_HELP_MESSAGE = """## Sequon Language

**Valid Residues:** 
- ARNDCEQGHILKMFPSTWYV
- X denoates all residues

**Special charaters:**
- "|" - deliminator between residue index
    - pre: -n|...|-3|-2|-1
    - post: +1|+2|+3|..|+n
- "!" - inverse operator
    - X!ARN - all residues **BUT** A, R, and N

## Example:

post: **X!P|ST** 
Allow modification if the next residue is any amino acid but P, and the next next residue is either S or T
- PEAT(HexNAc)LSIDE -> Allowed
- PERS(HexNAc)PSIDE -> Unallowed (P not permitted)

pre: **ARN|X!WYS**
Allow modification if the prior residue is anything but W, Y, and Z  and the prior prior residue is A,R, or N
- PEAT(HexNAc)LSIDE -> Allowed
- PERS(HexNAc)PSIDE -> Unallowed (S not permitted)"""
STATIC_MOD_HELP_COLUMN = "Static Modifications to use in search"
VARIABLE_MOD_HELP_COLUMN = "Variable Modifications to use in search"
MAX_VARIABLE_MODIFICATIONS_HELP_MESSAGE = "Number of variable modifications to allow per peptide"