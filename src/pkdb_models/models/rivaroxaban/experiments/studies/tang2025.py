"""Simulation experiment of the study Tang2025."""

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


class Tang2025(RivaroxabanSimulationExperiment):
    """Simulation experiment of Tang2025.

    Bioequivalence of a test (T) and a reference (R) tablet of rivaroxaban
    20 mg in two cohorts of 36 healthy Chinese subjects, fasted and fed (high-fat
    high-calorie meal 30 min before the dose), in a fully replicated four period
    crossover (TRTR/RTRT): every subject took each tablet twice (T1, R1, T2,
    R2), every curve is the mean of the same 36 subjects (Fig. 1, means without
    dispersion). The times of the copy are the nominal sampling times; Fed T1
    was read from the 0-12 hr inset (the 3 hr point and the tail after 12 hr
    are hidden by the other curves).
    """

    formulations: ClassVar[list[str]] = ["R1", "T1", "R2", "T2"]
    formulation_labels: ClassVar[dict[str, str]] = {
        "R1": "reference 20 mg (1st)",
        "T1": "test 20 mg (1st)",
        "R2": "reference 20 mg (2nd)",
        "T2": "test 20 mg (2nd)",
    }
    markers: ClassVar[dict[str, str]] = {"R1": "s", "T1": "^", "R2": "D", "T2": "v"}
    fasting_states: ClassVar[list[str]] = ["Fasting", "Fed"]
    fasting_labels: ClassVar[dict[str, str]] = {
        "Fasting": "fasted",
        "Fed": "fed",
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "Fasting": 64.2,  # kg
        "Fed": 65.6,  # kg
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        # the test and the reference tablet are simulated alike
        simulations = {}
        for fasting_state in self.fasting_states:
            fed = self.fasting_map[self.fasting_labels[fasting_state]]
            # the last samples are at 48 hr
            simulations[self.to_sid(f"po_{fasting_state}")] = Simulation(
                time_unit="min",
                end=50 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[fasting_state], "kg"),
                    "PODOSE_riv": Q(20, "mg"),
                    "GU__fed": Q(fed, "dimensionless"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for fasting_state in self.fasting_states:
            for formulation in self.formulations:
                key = f"{fasting_state}_{formulation}"
                mappings[self.to_sid(f"fm_po_{key}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"rivaroxaban_{key}"),
                        xid="time",
                        yid="mean",
                        yid_sd=None,
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

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for kp, fasting_state in enumerate(self.fasting_states):
            state = self.fasting_labels[fasting_state]
            color = self.fasting_colors[state]
            plots[kp].set_title(state)
            plots[kp].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
            # simulation
            plots[kp].add_data(
                task=self.to_sid(f"task_po_{fasting_state}"),
                xid="time",
                yid="[Cve_riv]",
                label=f"20 mg ({state})",
                color=color,
            )
            # data
            for formulation in self.formulations:
                plots[kp].add_data(
                    dataset=self.to_sid(f"rivaroxaban_{fasting_state}_{formulation}"),
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=self.formulation_labels[formulation],
                    color=color,
                    marker=self.markers[formulation],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Tang2025, output_dir=Tang2025.__name__)
