from pathlib import Path

from pkdb_models.models.data import collect_tsv_files

def collect_rivaroxaban_data():
    common_parent: Path = Path(__file__).parents[5]
    source_dir = common_parent / "pkdb_data" / "studies" / "rivaroxaban"
    target_dir = Path(__file__).parent / "rivaroxaban"

    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)

    # collect  (apixaban)
    def is_apixaban(study_name) -> bool:
        return study_name in ["Kreutz2017", "Frost2014"]

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "apixaban",
        target_dir=Path(__file__).parent / "apixaban",
        filter_study=is_apixaban,
    )

    # collect  (edoxaban)
    def is_edoxaban(study_name) -> bool:
        return study_name in ["Rohr2024", "Lenard2024", "Lenard2025"]

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "edoxaban",
        target_dir=Path(__file__).parent / "edoxaban",
        filter_study=is_edoxaban,
    )

if __name__ == "__main__":
    collect_rivaroxaban_data()

