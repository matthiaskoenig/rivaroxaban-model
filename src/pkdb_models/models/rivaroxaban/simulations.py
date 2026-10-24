"""The simulation experiments of the rivaroxaban model and their groups.

They are the registry of `MODEL`, the experiments its commands simulate.
"""

from sbmlsim.experiment import SimulationExperiment

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.misc import (
    CoagulationExperiment,
    DoseDependencyExperiment,
)
from pkdb_models.models.rivaroxaban.experiments.scans.scan_parameters import (
    RivaroxabanParameterScan,
)
from pkdb_models.models.rivaroxaban.experiments.studies import (
    Alalawneh2024,
    Ding2020,
    Frost2014,
    GenisNajera2022,
    Graff2007,
    Greenblatt2018,
    Jiang2010,
    Kreutz2017,
    Kubitza2005,
    Kubitza2005a,
    Kubitza2006,
    Kubitza2013,
    Lenard2024,
    Lenard2025,
    Liu2024,
    Metze2021,
    Moore2014,
    Rohr2024,
    Santos2024,
    Tang2025,
    Weinz2009,
    Zhao2009,
)

EXPERIMENTS: dict[str, list[type[SimulationExperiment]]] = {
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
        Metze2021,  # renal impairment
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
        Metze2021,
        Moore2014,
    ],
    "misc": [
        CoagulationExperiment,
        DoseDependencyExperiment,
    ],
    "scans": [
        RivaroxabanParameterScan,
    ],
}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scans"]


if __name__ == "__main__":
    MODEL.run_simulation_experiments(EXPERIMENTS["all"], MODEL.simulation_dir / "all")
