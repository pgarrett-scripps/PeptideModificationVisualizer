from pathlib import Path

import streamlit as st
from senpy.sqt.parser import read_file

st.title("SQT Editor")
st.write("Sets all mobility / timscore values to zero")


sqt_file = st.file_uploader("SQT File", ".sqt")

if sqt_file:

    with st.spinner("Reading SQT File"):
        sqt_lines = sqt_file.getvalue().decode('utf-8').split('\n')
        h_lines, s_lines = read_file(sqt_lines)

    for s_line in s_lines:
        for m_line in s_line.m_lines:
            m_line.tims_score = 0
            m_line.predicted_ook0 = 0

    with st.spinner("Writing SQT File"):
        lines = []
        for h_line in h_lines:
            lines.append(h_line.serialize(version='v2.1.0_ext'))
        for s_line in s_lines:
            lines.append(s_line.serialize(version='v2.1.0_ext'))
            for m_line in s_line.m_lines:
                lines.append(m_line.serialize(version='v2.1.0_ext'))
                for l_line in m_line.l_lines:
                    lines.append(l_line.serialize(version='v2.1.0_ext'))
        lines.append('\n')

    st.download_button("Download Zero SQT", "".join(lines), file_name=f"{Path(sqt_file.name).stem}_ZERO.sqt")

