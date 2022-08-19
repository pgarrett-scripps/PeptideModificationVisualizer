from io import StringIO

from pyopenms import FASTAEntry


def convert_fasta(entries):
    lines = []
    for fasta_entry in entries:
        lines.append(f">{fasta_entry.identifier}")
        lines.append(f"{fasta_entry.sequence}")
    return "\n".join(lines)


def load_fasta(file):
    stringio = StringIO(file.getvalue().decode("utf-8"))
    entries = []
    for line in stringio:
        if line[0] == ">":  # new protein
            fasta_entry = FASTAEntry()
            entries.append(fasta_entry)

            identifier = line.rstrip().split(" ")[0][1:]
            entries[-1].identifier = identifier
        else:  # protein sequence
            entries[-1].sequence += line.rstrip()
    return entries
