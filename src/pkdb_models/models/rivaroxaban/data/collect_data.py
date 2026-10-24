"""Collect the data of the rivaroxaban studies from the PK-DB study folders.

The studies are copied from a `pkdb_data` checkout (`PKDB_DATA`, else
`../pkdb_data`, see `pkpdsim.data.studies_dir`) into the data folder of the
model, see `pkpdsim.data.collect_tsv_files`.
Besides the rivaroxaban studies, studies of other factor Xa inhibitors have
rivaroxaban groups: Frost2014 and Kreutz2017 (apixaban), Lenard2024,
Lenard2025, Metze2021 and Rohr2024 (edoxaban).
"""

from pathlib import Path

from pkpdsim.data import collect_tsv_files
from pkpdsim.data import studies_dir as pkdb_studies_dir

from pkdb_models.models.rivaroxaban import (
    DATA_PATH_APIXABAN,
    DATA_PATH_EDOXABAN,
    DATA_PATH_RIVAROXABAN,
)

#: the studies of other substances which the rivaroxaban experiments use
APIXABAN_STUDIES = {"Frost2014", "Kreutz2017"}
EDOXABAN_STUDIES = {"Lenard2024", "Lenard2025", "Metze2021", "Rohr2024"}


def collect_data(studies_dir: Path | None = None) -> None:
    """Copy the data of the rivaroxaban studies into the data folder of the model.

    Args:
        studies_dir: the `studies` folder of pkdb_data, by default the one of
            the checkout `pkpdsim.data.studies_dir` finds.
    """
    if studies_dir is None:
        studies_dir = pkdb_studies_dir()
    collect_tsv_files(
        source_dir=studies_dir / "rivaroxaban", target_dir=DATA_PATH_RIVAROXABAN
    )
    collect_tsv_files(
        source_dir=studies_dir / "apixaban",
        target_dir=DATA_PATH_APIXABAN,
        filter_study=lambda study: study in APIXABAN_STUDIES,
    )
    collect_tsv_files(
        source_dir=studies_dir / "edoxaban",
        target_dir=DATA_PATH_EDOXABAN,
        filter_study=lambda study: study in EDOXABAN_STUDIES,
    )


if __name__ == "__main__":
    collect_data()
