"""Run all simulation experiments."""
import shutil
from pathlib import Path
from typing import List

from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.helpers import run_experiments
from pkdb_models.models.rivaroxaban.experiments.studies import *
from pkdb_models.models.rivaroxaban.experiments.misc import *
from pkdb_models.models.rivaroxaban.experiments.scans.scan_parameters import RivaroxabanParameterScan
from sbmlutils import log
from sbmlsim.plot import Figure

Figure.legend_fontsize = 8
Figure.fig_dpi = 300


logger = log.get_logger(__name__)

EXPERIMENTS = {
    "studies": [
        Alalawneh2024,
        Ding2020,  # food
        Frost2014,
        GenisNajera2022,
        Graff2007,  # dose-dependency
        Greenblatt2018,
        Jiang2010,  # dose-dependency
        Kreutz2017,  # food
        Kubitza2005,  # dose-dependency
        Kubitza2005a,  # dose-dependency, application-form
        Kubitza2006,  # food
        Kubitza2013,
        Lenard2024,
        Lenard2025,
        Liu2024,  # food
        Moore2014,
        Rohr2024,
        Santos2024,
        Tang2025,
        Weinz2009,
        Zhao2009,  # dose-dependency
    ],
    "pharmacodynamics": [
        Alalawneh2024,
        Graff2007,
        Jiang2010,
        Kreutz2017,
        Kubitza2005,
        Kubitza2005a,
        Zhao2009,
    ],
    "bodyweight": [
        Alalawneh2024,
    ],
    "dose_dependency": [
        Graff2007,
        Jiang2010,
        Kubitza2005,
        Kubitza2005a,
        Zhao2009,
    ],
    "food": [
        Kreutz2017,
        Ding2020,
        Kubitza2006,
        Liu2024,
        Santos2024,
    ],
    "hepatic_impairment": [
        Kubitza2013,
    ],
    "renal_impairment": [
        Greenblatt2018,
        Moore2014,
    ],
    "misc": [
        DoseDependencyExperiment
    ],
    "scan": [
        RivaroxabanParameterScan,
    ]

}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["scan"] + EXPERIMENTS["misc"]

def run_simulation_experiments(
    selected: str = None,
    experiment_classes: List = None,
    output_dir: Path = None
) -> None:
    """Run rivaroxaban simulation experiments."""

    Figure.fig_dpi = 600
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = rivaroxaban.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = rivaroxaban.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")



if __name__ == "__main__":
    """Run experiments."""

    # selected = "all"
    # selected = "studies"
    # selected = "pharmacodynamics"

    run_simulation_experiments(selected="studies")
