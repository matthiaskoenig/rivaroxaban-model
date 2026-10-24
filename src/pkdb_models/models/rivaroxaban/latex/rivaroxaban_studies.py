"""LaTeX table of the rivaroxaban studies.

Writes `rivaroxaban_studies.tex` from the overview of the studies in
`rivaroxaban_studies.tsv`.
"""

from pathlib import Path

import pandas as pd

#: overview of the studies
TSV_PATH = Path(__file__).parent / "rivaroxaban_studies.tsv"
#: the generated table
LATEX_PATH = Path(__file__).parent / "rivaroxaban_studies.tex"

#: preamble of the table, up to the column headers; some lines end with a
#: space, which is kept as it is part of the generated table
LATEX_HEADER = "\n".join(
    [
        r"%\begin{landscape}",
        r"\begin{table}[ht!]",
        r"\centering",
        r"\label{tab:rivaroxaban_studies}",
        r"\tabcolsep=3.5pt\relax",
        r"\scriptsize",
        r"\begin{threeparttable} ",
        "",
        r"\caption{\textbf{Overview of rivaroxaban studies.}  ",
        r"This table summarizes the included studies with key identifiers (\emph{Study}, \emph{PMID}, \emph{PK-DB ID}), dosing regimen (\emph{Dosing}), dose (\emph{Dose}), and participant characteristics: \emph{H} (healthy), \emph{RI} (renal impairment), and \emph{HI} (hepatic impairment). Fasting status is indicated as \emph{FA} (fasted) or \emph{FE} (fed). Pharmacokinetic data include \emph{RP} (plasma rivaroxaban), \emph{RU} (urine rivaroxaban), and \emph{RF} (feces rivaroxaban); pharmacodynamic data include \emph{Xa} (Factor Xa inhibition), \emph{PT} (prothrombin time), and \emph{aP} (activated partial thromboplastin time).}",
        "",
        "",
        r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}} ",
        r"  p{1.8cm}  % Study",
        r"  p{1.3cm}  % PK-DB",
        r"  p{0.9cm}  % Dosing",
        r"  p{1.2cm}  % Dose [mg]",
        r"  p{0.3cm}  % H",
        r"  p{0.3cm}  % RI",
        r"  p{0.3cm}  % HI",
        r"  p{0.3cm}  % fasted",
        r"  p{0.3cm}  % fed",
        r"  p{0.3cm}  % riv plasma",
        r"  p{0.3cm}  % riv urine",
        r"  p{0.3cm}  % riv feces",
        r"  p{0.3cm}  % Xa",
        r"  p{0.3cm}  % PT",
        r"  p{0.3cm}  % aPTT",
        r"  p{0.4cm}  % ref",
        r"}",
        r"\toprule",
    ]
)


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""
    # Selecting needed columns
    df = df.loc[
        :,
        [
            "study",
            "pkdb",
            "dosing",
            "dose",
            "healthy",
            "renal impairment",
            "hepatic impairment",
            "fasted",
            "fed",
            "riv plasma",
            "riv urine",
            "riv feces",
            "Xa",
            "PT",
            "aPTT",
            "pmid",
        ],
    ].copy()

    # Renaming columns
    df = df.rename(
        columns={
            "study": "Study",
            "pkdb": "PK-DB",
            "dosing": "Dosing",
            "dose": "Dose [mg]",
            "healthy": "H",
            "renal impairment": "RI",
            "hepatic impairment": "HI",
            "fasted": "FA",
            "fed": "FE",
            "riv plasma": "RP",
            "riv urine": "RU",
            "riv feces": "RF",
            "aPTT": "aP",
            "pmid": "S",
        }
    )
    # Remove pipe characters from "Dose [mg]" column
    if "Dose [mg]" in df.columns:
        df["Dose [mg]"] = df["Dose [mg]"].str.replace("|", ",", regex=False)

    # Convert True/False to checkmarks
    checkmark_columns = [
        "H",
        "RI",
        "HI",
        "FA",
        "FE",
        "RP",
        "RU",
        "RF",
        "Xa",
        "PT",
        "aP",
    ]

    for col in checkmark_columns:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: r"\checkmark" if str(x).strip().upper() == "TRUE" else ""
            )

    # Column headers row
    column_names = df.columns.tolist()
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in column_names]) + r" \\"

    # LaTeX body
    latex_body = f"{column_headers}\n\\midrule\n"

    for _, row in df.iterrows():
        values = list(row.astype(str).values)

        # references
        study_col_index = 0
        ref_col_index = 15
        values[ref_col_index] = f"\\cite{{{values[study_col_index]}}}"

        # pkdb urls
        pkdb_col_index = 1
        values[pkdb_col_index] = (
            f"\\href{{https://identifiers.org/pkdb:{values[pkdb_col_index]}}}"
            f"{{{values[pkdb_col_index]}}}"
            if values[pkdb_col_index]
            else ""
        )

        row_str = " & ".join(v.replace("TRUE", r"\checkmark") for v in values) + r" \\"

        latex_body += row_str + "\n"

    latex_footer = r"""
\bottomrule
\end{tabular*}

\end{threeparttable}
\end{table}
%\end{landscape}
""".strip("\n")

    # Combine all parts
    return LATEX_HEADER + "\n" + latex_body + latex_footer


def main() -> None:
    """Write the LaTeX table of the studies."""
    # Load TSV file
    df = pd.read_csv(TSV_PATH, sep="\t", header=1, nrows=14)

    # Replace NaN values with an empty string
    df = df.fillna("")

    # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with LATEX_PATH.open("w") as f_tex:
        f_tex.write(latex_str)


if __name__ == "__main__":
    main()
