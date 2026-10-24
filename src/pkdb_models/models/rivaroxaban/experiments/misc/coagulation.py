"""Effect of rivaroxaban on the coagulation model alone."""

from typing import ClassVar, override

import numpy as np
from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.coagulation_experiment import (
    CoagulationSimulationExperiment,
)


class CoagulationExperiment(CoagulationSimulationExperiment):
    """Effect of rivaroxaban."""

    riv_values: ClassVar[np.ndarray] = np.linspace(
        0.0, 300.0 / 435.88, num=6
    )  # [µmole/l]  # FIXME values not right
    colors: ClassVar[list[str]] = [
        "black",
        "#c6dbef",
        "#9ecae1",
        "#6baed6",
        "#3182bd",
        "#08519c",
    ]

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for k, riv in enumerate(self.riv_values):
            # rivaroxaban for 24 hr after a baseline of 5 hr, then 24 hr without
            simulations[f"coagulation_{k}"] = Simulation(
                time_unit="min",
                end=(5 + 24 + 24) * 60,
                preinit_changes={
                    **self.default_changes(),
                    "[Cve_riv]": Q(0, "µM"),
                },
                changes=[
                    Change(5 * 60, {"[Cve_riv]": Q(riv, "µM")}),
                    Change((5 + 24) * 60, {"[Cve_riv]": Q(0, "µM")}),
                ],
                time_shift=-5 * 60,
            )

        return simulations

    @override
    def figures(self) -> dict[str, Figure]:
        # the concentration of rivaroxaban and its effects, by panel
        info = [
            ("[Cve_riv]", 0),
            ("Xa_inhibition", 1),
            ("PT", 2),
            ("aPTT", 3),
            ("PT_change", 4),
            ("aPTT_change", 5),
            ("PT_ratio", 6),
            ("aPTT_ratio", 7),
        ]

        fig_id = "Fig_coagulation"
        fig = Figure(
            experiment=self,
            sid=fig_id,
            num_rows=4,
            num_cols=2,
            name="Coagulation",
        )
        plots = fig.create_plots(xaxis=Axis("time", unit="hour"), legend=True)

        for sid, ksid in info:
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for k, riv in enumerate(self.riv_values):
                plots[ksid].add_data(
                    task=f"task_coagulation_{k}",
                    xid="time",
                    yid=sid,
                    label=f"{riv:.2f} µM",
                    color=self.colors[k],
                )

        return {fig_id: fig}


if __name__ == "__main__":
    MODEL.run_experiments(
        CoagulationExperiment, output_dir=CoagulationExperiment.__name__
    )
