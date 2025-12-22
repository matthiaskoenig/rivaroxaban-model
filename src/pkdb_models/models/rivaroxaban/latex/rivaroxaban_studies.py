from pathlib import Path
import pandas as pd


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    # Preamble:
    latex_header = r"""
%\begin{landscape}
\begin{table}[ht!]
\centering
\label{tab:rivaroxaban_studies}
\tabcolsep=3.5pt\relax
\scriptsize
\begin{threeparttable} 

\caption{\textbf{Overview of rivaroxaban studies.}  
This table summarizes the included studies with key identifiers (\emph{Study}, \emph{PMID}, \emph{PK-DB ID}), dosing regimen (\emph{Dosing}), dose (\emph{Dose}), and participant characteristics: \emph{H} (healthy), \emph{RI} (renal impairment), and \emph{HI} (hepatic impairment). Fasting status is indicated as \emph{FA} (fasted) or \emph{FE} (fed). Pharmacokinetic data include \emph{RP} (plasma rivaroxaban), \emph{RU} (urine rivaroxaban), and \emph{RF} (feces rivaroxaban); pharmacodynamic data include \emph{Xa} (Factor Xa inhibition), \emph{PT} (prothrombin time), and \emph{aP} (activated partial thromboplastin time).}


\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}} 
  p{1.8cm}  % Study
  p{1.3cm}  % PK-DB
  p{0.9cm}  % Dosing
  p{1.2cm}  % Dose [mg]
  p{0.3cm}  % H
  p{0.3cm}  % RI
  p{0.3cm}  % HI
  p{0.3cm}  % fasted
  p{0.3cm}  % fed
  p{0.3cm}  % riv plasma
  p{0.3cm}  % riv urine
  p{0.3cm}  % riv feces
  p{0.3cm}  % Xa
  p{0.3cm}  % PT
  p{0.3cm}  % aPTT
  p{0.4cm}  % ref
}
\toprule
""".strip('\n')

    # Selecting needed columns
    df = df.loc[:, [
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
                   ]].copy()

    # Renaming columns
    df = df.rename(columns={
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
    })
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
        values[ref_col_index] = (
            f"\\cite{{{values[study_col_index]}}}"
        )

        # # pubmed urls
        # pmid_col_index = 1
        # values[pmid_col_index] = (
        #     f"\\href{{https://pubmed.ncbi.nlm.nih.gov/{values[pmid_col_index]}}}"
        #     f"{{{values[pmid_col_index]}}}"
        #     if values[pmid_col_index] else ""
        # )

        # pkdb urls
        pkdb_col_index = 1
        values[pkdb_col_index] = (
            f"\\href{{https://identifiers.org/pkdb:{values[pkdb_col_index]}}}"
            f"{{{values[pkdb_col_index]}}}"
            if values[pkdb_col_index] else ""
        )

        row_str = " & ".join(
            v.replace("TRUE", r"\checkmark")
            for v in values
        ) + r" \\"

        latex_body += row_str + "\n"

    latex_footer = r"""
\bottomrule
\end{tabular*}

\end{threeparttable}
\end{table}
%\end{landscape}
""".strip('\n')

    # Combine all parts
    full_latex = latex_header + "\n" + latex_body + latex_footer
    return full_latex


if __name__ == "__main__":
    # Input and output file paths
    tsv_path = Path(__file__).parent / 'rivaroxaban_studies.tsv'
    latex_path = Path(__file__).parent / 'rivaroxaban_studies.tex'

    # Load TSV file
    df = pd.read_csv(tsv_path, sep="\t", header=1, nrows=14)

    # Replace NaN values with an empty string
    df = df.fillna("")

    # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)
