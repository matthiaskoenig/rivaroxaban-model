"""Simulation experiment of the study Ding2020."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)


class Ding2020(RivaroxabanSimulationExperiment):
    """Simulation experiment of Ding2020.

    Bioequivalence of a generic (T) and Xarelto (R) 10 mg tablet in healthy
    Chinese volunteers, a fasted part (overnight fast of 10 hr) and a fed part
    (high-fat, high-calorie breakfast), each a replicate crossover of 28
    subjects (27 completed the fed part). T and R are the same subjects and
    are simulated alike; the curves are means +- SD (Fig1 fasted, Fig2 fed).
    """

    formulations: ClassVar[list[str]] = ["T10", "R10"]
    fasting_states: ClassVar[list[str]] = ["Fasting", "Fed"]
    bodyweights: ClassVar[dict[str, float]] = {
        "Fasting": 63.0,  # kg
        "Fed": 62.0,  # kg
    }
    fed: ClassVar[dict[str, float]] = {
        "Fasting": RivaroxabanSimulationExperiment.fasting_map["fasted"],
        "Fed": RivaroxabanSimulationExperiment.fasting_map["fed"],
    }
    colors: ClassVar[dict[str, str]] = {
        "Fasting": RivaroxabanSimulationExperiment.fasting_colors["fasted"],
        "Fed": RivaroxabanSimulationExperiment.fasting_colors["fed"],
    }
    fasting_labels: ClassVar[dict[str, str]] = {
        "Fasting": "fasted",
        "Fed": "fed",
    }
    formulation_labels: ClassVar[dict[str, str]] = {
        "T10": "test 10 mg",
        "R10": "Xarelto 10 mg",
    }
    markers: ClassVar[dict[str, str]] = {
        "T10": "^",
        "R10": "o",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for fasting_state in self.fasting_states:
            tcsims[self.to_sid(f"po_{fasting_state}")] = Simulation(
                time_unit="min",
                end=50 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[fasting_state], "kg"),
                    "PODOSE_riv": Q(10, "mg"),
                    "GU__fed": Q(self.fed[fasting_state], "dimensionless"),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for fasting_state in self.fasting_states:
            for formulation in self.formulations:
                mappings[self.to_sid(f"fm_po_{fasting_state}_{formulation}")] = (
                    FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(
                                f"rivaroxaban_{fasting_state}_{formulation}"
                            ),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=self.to_sid(f"task_po_{fasting_state}"),
                            xid="time",
                            yid="[Cve_riv]",
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.FASTED
                            if fasting_state == "Fasting"
                            else Fasting.FED,
                            coadministration=Coadministration.NONE,
                        ),
                    )
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1_2",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy), 10 mg po",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=2e-4
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False

        for plot in plots:
            for fasting_state in self.fasting_states:
                fasting_label = self.fasting_labels[fasting_state]
                # simulation
                plot.add_data(
                    task=self.to_sid(f"task_po_{fasting_state}"),
                    xid="time",
                    yid="[Cve_riv]",
                    label=fasting_label,
                    color=self.colors[fasting_state],
                )

                # data
                for formulation in self.formulations:
                    plot.add_data(
                        dataset=self.to_sid(
                            f"rivaroxaban_{fasting_state}_{formulation}"
                        ),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=(
                            f"{fasting_label}, {self.formulation_labels[formulation]}"
                        ),
                        color=self.colors[fasting_state],
                        marker=self.markers[formulation],
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Ding2020, output_dir=Ding2020.__name__)
