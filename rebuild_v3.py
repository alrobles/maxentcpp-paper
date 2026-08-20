#!/usr/bin/env python3
"""Rebuild paper_v3.tex from paper_v3.md using repo's existing paper_v3.tex template."""
import re, subprocess, os

REPO = "/home/alrobles/maxentcpp-paper"
MD = f"{REPO}/paper_v3.md"
PANDOC_TEX = "/tmp/paper_v3_pandoc.tex"
OUT = f"{REPO}/paper_v3.tex"
CURRENT = f"{REPO}/paper_v3.tex"
INCLUDE_APPENDIX = True

# 1. pandoc conversion
subprocess.run(["pandoc", MD, "-o", PANDOC_TEX, "--natbib", "--bibliography", f"{REPO}/paper.bib"],
               check=True)
with open(PANDOC_TEX) as fh:
    tex = fh.read()
with open(CURRENT) as fh:
    cur = fh.read()

header_end = cur.index("\\section*{Summary}")
header = cur[:header_end]

body_start = tex.index("\\section{Summary}")
body = tex[body_start:]

body = re.sub(r"\\citep\{([^}]*)\}", r"\\cite{\1}", body)
body = re.sub(r"\\citet\{([^}]*)\}", r"\\cite{\1}", body)
body = re.sub(r"\\cite\{([^}]*)\}",
              lambda m: "\\cite{" + m.group(1).replace(" ", "") + "}", body)
body = re.sub(r"\\section\{([^}]*)\}\\label\{[^}]*\}", r"\\section*{\1}", body)
body = re.sub(r"\\subsection\{([^}]*)\}\\label\{[^}]*\}", r"\\subsection*{\1}", body)
body = re.sub(r"\\subsubsection\{([^}]*)\}\\label\{[^}]*\}", r"\\subsubsection*{\1}", body)

out = header + body
if "\\end{document}" not in out:
    out += "\n\\nolinenumbers\n\n% bibliography\n\\bibliography{paper}\n\\bibliographystyle{abbrv}\n\n"
    if INCLUDE_APPENDIX:
        out += "% appendix (shared body with appendix_benchmarks.tex)\n\\appendix\n\\input{appendix_body}\n\n"
    out += "\\end{document}\n"

# The appendix_body is NOT part of the pandoc body; keep the template's appendix close.
# If lambda-file section header reference issues, handled by compile.

with open(OUT, "w") as fh:
    fh.write(out)
print(f"wrote {OUT} ({len(out)} chars)")
print("citep remaining:", out.count("\\citep"))
print("citet remaining:", out.count("\\citet"))
print("numbered sections:", len(re.findall(r"\\(?:sub)*section\{", out)))
print("undef labels:", len(re.findall(r"\\label\{[^}]*\}", out)))