"""Dependency of the pharmacokinetics and pharmacodynamics on the dose."""

from typing import ClassVar, override

from sbmlsim import Q
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Change, Simulation

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)


class DoseDependencyExperiment(RivaroxabanSimulationExperiment):
    """Single and multiple doses of rivaroxaban (IV, PO)."""

    doses: ClassVar[list[int]] = [0, 10, 20, 40, 80]  # [mg]
    colors: ClassVar[list[str]] = [
        "black",
        "tab:orange",
        "tab:blue",
        "tab:red",
        "tab:green",
    ]
    routes: ClassVar[list[str]] = ["IV", "PO"]
    dosing: ClassVar[list[str]] = ["single", "multiple"]

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}

        for route in self.routes:
            for dose in self.doses:
                for dosing in self.dosing:
                    # single dosing
                    if dosing == "single":
                        simulations[f"riv_{route}_{dose}_{dosing}"] = Simulation(
                            time_unit="min",
                            end=24 * 60,
                            preinit_changes={
                                **self.default_changes(),
                                f"{route}DOSE_riv": Q(dose, "mg"),
                            },
                        )
                    elif dosing == "multiple":
                        # a dose every 24 hr for 7 days, the urine is
                        # collected from dose to dose; a dose adds to what is
                        # left of the earlier ones
                        times = [24 * 60 * k for k in range(1, 7)]
                        simulations[f"riv_{route}_{dose}_{dosing}"] = Simulation(
                            time_unit="min",
                            end=7 * 24 * 60,
                            preinit_changes={
                                **self.default_changes(),
                                f"{route}DOSE_riv": Q(dose, "mg"),
                            },
                            changes=[
                                Change(
                                    times,
                                    {
                                        **self.add_doses(
                                            times, {f"{route}DOSE_riv": Q(dose, "mg")}
                                        ).values,
                                        # reset urinary amounts
                                        "Aurine_riv": Q(0, "mmole"),
                                        "Aurine_rx": Q(0, "mmole"),
                                    },
                                )
                            ],
                        )

        return simulations

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> dict[str, Figure]:
        """Timecourses of the pharmacokinetics and the effects by route and dosing."""
        figures = []

        for route in self.routes:
            for dosing in self.dosing:
                fig = Figure(
                    experiment=self,
                    sid=f"Fig_dose_dependency_pk_{route}_{dosing}",
                    num_rows=6,
                    num_cols=3,
                    name=f"Dose dependency rivaroxaban ({route}, {dosing})",
                    settings={"legend_scope": "row"},
                )

                plots = fig.create_plots(
                    xaxis=Axis("time", unit="hr" if dosing == "single" else "day"),
                    legend=True,
                )
                infos = [
                    # plasma
                    ("[Cve_riv]", 0),
                    ("[Cve_rx]", 1),
                    ("[Cve_riv_rx]", 2),
                    # urine
                    ("Aurine_riv", 3),
                    ("Aurine_rx", 4),
                    ("Aurine_riv_rx", 5),
                    # feces
                    ("Afeces_riv", 6),
                    ("Afeces_rx", 7),
                    ("Afeces_riv_rx", 8),
                    ("PT", 9),
                    ("PT_change", 10),
                    ("PT_ratio", 11),
                    ("aPTT", 12),
                    ("aPTT_change", 13),
                    ("aPTT_ratio", 14),
                    ("Xa_inhibition", 15),
                    ("KI__RIVEX", 16),
                    (f"{route}DOSE_riv", 17),
                ]
                for sid, ksid in infos:
                    plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

                for sid, ksid in infos:
                    for kval, dose in enumerate(self.doses):
                        plots[ksid].add_data(
                            task=f"task_riv_{route}_{dose}_{dosing}",
                            xid="time",
                            yid=sid,
                            label=f"{dose} mg",
                            color=self.colors[kval],
                        )

                figures.append(fig)
        return self.figures_by_sid(*figures)


if __name__ == "__main__":
    MODEL.run_experiments(
        DoseDependencyExperiment, output_dir=DoseDependencyExperiment.__name__
    )
