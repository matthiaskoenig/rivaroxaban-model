from typing import Dict

import numpy as np
from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rivaroxaban.experiments.coagulation_experiment import CoagulationSimulationExperiment
from pkdb_models.models.rivaroxaban.helpers import run_experiments


class CoagulationExperiment(CoagulationSimulationExperiment):
    """Effect of rivaroxaban."""

    riv_values = np.linspace(0.0, 300.0/435.88, num=6)  # [µmole/l]  # FIXME values not right
    colors = ['black', '#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c']

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        baseline_changes = {}

        for k, riv in enumerate(self.riv_values):
            tcsims[f"coagulation_{k}"] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=5 * 60,  # [min] # simulate 1 day
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "[Cve_riv]": Q_(0, "µM"),
                        **baseline_changes,
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[Cve_riv]": Q_(riv, "µM")
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[Cve_riv]": Q_(0, "µM")
                    },
                ),
                ],
                time_offset=-5*60
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        info = [
            ("PT", 0),
            ("PT_change", 1),
            ("PT_ratio", 2),

            ("aPTT", 3),
            ("aPTT_change", 4),
            ("aPTT_ratio", 5),

            ("Xa_inhibition", 6),
        ]

        fig = Figure(
            experiment=self,
            sid=f"Fig_coagulation",
            num_rows=3,
            num_cols=3,
            name=f"Coagulation",
        )
        plots = fig.create_plots(xaxis=Axis("time", unit="hour"), legend=True)

        for sid, ksid in info:
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for k, riv in enumerate(self.riv_values):
                plots[ksid].add_data(
                    task=f"task_coagulation_{k}",
                    xid="time",
                    yid=sid,
                    label=f"{riv:.1f} nM",
                    color=self.colors[k],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(
       CoagulationExperiment, output_dir=CoagulationExperiment.__name__
    )
