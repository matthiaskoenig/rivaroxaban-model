"""Rivaroxaban: PBPK/PD model of rivaroxaban and its anticoagulant effect.

The paths of the model and its model object `MODEL`, which the commands
`run_simulation`, `run_fit`, `run_report` and `run_identifiability` of pkpdsim
run by its package `pkdb_models.models.rivaroxaban`.
"""

from pathlib import Path

from pkpdsim.model import FigureSettings, PKPDModel, SimulatorSettings
from pkpdsim.release.settings import Author, ReleaseSettings

from pkdb_models.models.rivaroxaban.tolerances import ABSOLUTE_TOLERANCE

RIVAROXABAN_PATH = Path(__file__).parent

MODEL_BASE_PATH = RIVAROXABAN_PATH / "models" / "output" / "models"
MODEL_PATH = MODEL_BASE_PATH / "rivaroxaban_body_flat.xml"

DATA_PATH_BASE = RIVAROXABAN_PATH / "data"
DATA_PATH_RIVAROXABAN = DATA_PATH_BASE / "rivaroxaban"
DATA_PATH_APIXABAN = DATA_PATH_BASE / "apixaban"
DATA_PATH_EDOXABAN = DATA_PATH_BASE / "edoxaban"
DATA_PATHS = [
    DATA_PATH_RIVAROXABAN,
    DATA_PATH_APIXABAN,
    DATA_PATH_EDOXABAN,
]

MODEL = PKPDModel(
    name="rivaroxaban",
    package=__name__,
    title="Rivaroxaban PBPK/PD model",
    path=RIVAROXABAN_PATH,
    model_path=MODEL_PATH,
    data_paths=tuple(DATA_PATHS),
    archive_path=MODEL_BASE_PATH.parent / "rivaroxaban_model.omex",
    factory=".models.factory:create_rivaroxaban_models",
    experiments=".simulations:EXPERIMENTS",
    fit_definitions=".fitting.fitting:FIT_DEFINITIONS",
    # the adopted fits, the factory takes the fitted values from their records
    fits=("PK_20261018_075718__82026", "PD_20261018_082509__94e87"),
    simulator=SimulatorSettings(
        absolute_tolerance=ABSOLUTE_TOLERANCE, relative_tolerance=1e-10
    ),
    # the resolution of the PNG figures; at 600 dpi the large figures need up to
    # 12 GB of memory to render
    figures=FigureSettings(formats=("svg",), legend_fontsize=10, fig_dpi=300),
    fit_identifiability=True,
    # the release of the model repository
    release=ReleaseSettings(
        title=(
            "Physiologically based pharmacokinetic/pharmacodynamic (PBPK/PD) model"
            " of rivaroxaban"
        ),
        repository="matthiaskoenig/rivaroxaban-model",
        doi="10.5281/zenodo.15157735",
        authors=(
            Author(
                "Casabianca, Elisabetta",
                None,
                "Faculty of Pharmacy, Université Grenoble Alpes, Grenoble, France",
            ),
            Author(
                "Myshkina, Mariia",
                "0009-0001-1433-4293",
                "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver",
            ),
            Author(
                "König, Matthias",
                "0000-0003-1725-179X",
                (
                    "Humboldt-Universität zu Berlin, Faculty of Life Sciences, Department of Biology, Institute of Theoretical Biology, Systems Medicine of the Liver; University of Stuttgart, Institute of Structural Mechanics and Dynamics in Aerospace Engineering; University of Southern Denmark, Clinical Pharmacology, Pharmacy and Environmental Medicine, Department of Public Health; University Hospital Schleswig-Holstein, Campus Lübeck, First Department of Medicine, Lübeck, Germany"
                ),
            ),
        ),
        keywords=(
            "model",
            "rivaroxaban",
            "factor Xa inhibitor",
            "anticoagulation",
            "pharmacokinetics",
            "pharmacodynamics",
            "PBPK/PD",
            "physiologically based pharmacokinetic/ pharmacodynamic model",
            "SBML",
        ),
        docker_image="matthiaskoenig/rivaroxaban",
        funding=(
            "Matthias König was supported by the Federal Ministry of Research, Technology and Space (BMFTR, Germany) within ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B)."
        ),
        description=(
            "This repository provides the rivaroxaban physiologically based"
            " pharmacokinetic/pharmacodynamic (PBPK/PD) model: the dose and food"
            " dependent dissolution and the absorption, the tissue distribution,"
            " the hepatic metabolism (CYP3A4, CYP2J2, hydrolysis), the renal"
            " excretion by filtration and active secretion, the excretion of the"
            " metabolites and the anticoagulant effects (PT, aPTT, inhibition"
            " of factor Xa)."
        ),
        first_year=2025,
    ),
)
"""The rivaroxaban model: its commands, see `pkpdsim.model`."""
