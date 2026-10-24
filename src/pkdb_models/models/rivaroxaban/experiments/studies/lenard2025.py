"""Simulation experiment of the study Lenard2025."""

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
from pkdb_models.models.rivaroxaban.experiments.studies.lenard2024 import impute_sd


class Lenard2025(RivaroxabanSimulationExperiment):
    """Simulation experiment of Lenard2025.

    Microdose cocktail of apixaban 25 µg, rivaroxaban 25 µg and edoxaban 50 µg
    as an oral solution in three cohorts of healthy subjects (the trial of
    Lenard2024), each given the cocktail alone and at the steady state of a
    perpetrator: carbamazepine (CAR, Fig2, n=12), gabapentin (GAB, Fig3, n=11)
    and pregabalin (PRE, Fig4, n=12). Only the rivaroxaban data is used.
    """

    #: cohorts with the figure of their data
    cohorts: ClassVar[dict[str, str]] = {"car": "Fig2", "gab": "Fig3", "pre": "Fig4"}
    #: the intervention of the perpetrator arm of a cohort
    perpetrators: ClassVar[dict[str, str]] = {
        "car": "API25, RIV25, EDO50, CAR",
        "gab": "API25, RIV25, EDO50, GAB",
        "pre": "API25, RIV25, EDO50, PRE",
    }
    perpetrator_labels: ClassVar[dict[str, str]] = {
        "car": "carbamazepine",
        "gab": "gabapentin",
        "pre": "pregabalin",
    }
    coadministrations: ClassVar[dict[str, Coadministration]] = {
        "car": Coadministration.CARBAMAZEPINE,
        "gab": Coadministration.GABAPENTIN,
        "pre": Coadministration.PREGABALIN,
    }
    colors: ClassVar[dict[str, str]] = {
        "car": "tab:blue",
        "gab": "tab:orange",
        "pre": "tab:green",
    }
    #: the simulation of a perpetrator arm: gabapentin and pregabalin are
    #: neither metabolized by CYP nor inducers or inhibitors of CYP3A4 or P-gp
    #: (negative controls, observed AUC ratios 0.92 and 0.98), their arms are
    #: the simulation of the cocktail alone
    perpetrator_tasks: ClassVar[dict[str, str]] = {
        "car": "task_RIV25_CAR",
        "gab": "task_RIV25",
        "pre": "task_RIV25",
    }
    #: induction by carbamazepine at steady state (not fitted): a strong
    #: inducer of CYP3A4 and P-gp (PXR/CAR); the hepatic CYP3A4 activity 3 fold
    #: and the renal P-gp activity 2 fold give the AUC ratio of 0.54 for the
    #: observed 0.56 (Cmax 0.75) of this study (rifampicin with 20 mg: AUC
    #: 0.50, Mueck2013); without an intestinal CYP3A4 or P-gp the model does
    #: not describe the decrease of Cmax
    induction: ClassVar[dict[str, dict[str, float]]] = {
        "RIV25": {"LI__f_cyp3a4": 1.0, "KI__f_pgp": 1.0},
        "RIV25_CAR": {"LI__f_cyp3a4": 3.0, "KI__f_pgp": 2.0},
    }
    dose = 0.025  # [mg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for cohort, fig_id in self.cohorts.items():
            df = self.load_dataframe(fig_id)
            df = df[df.substance == "rivaroxaban"]
            for intervention, df_label in df.groupby("intervention"):
                key = (
                    f"riv25_{cohort}"
                    if intervention == self.perpetrators[cohort]
                    else f"riv25_{cohort}_baseline"
                )
                dset = DataSet.from_df(impute_sd(df_label), self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[key] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            key: Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "POSOLUTION_riv": Q(self.dose, "mg"),
                    **{
                        pid: Q(value, "dimensionless") for pid, value in changes.items()
                    },
                },
            )
            for key, changes in self.induction.items()
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for cohort in self.cohorts:
            for key, task, coadministration in [
                (f"riv25_{cohort}_baseline", "task_RIV25", Coadministration.NONE),
                (
                    f"riv25_{cohort}",
                    self.perpetrator_tasks[cohort],
                    self.coadministrations[cohort],
                ),
            ]:
                mappings[f"fm_{key}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=key,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(self, task=task, xid="time", yid="[Cve_riv]"),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=coadministration,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_Fig3_Fig4",
            num_cols=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        for kp, cohort in enumerate(self.cohorts):
            plots[kp].set_title(f"{self.perpetrator_labels[cohort]} cohort")
            plots[kp].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
            for key, task, label, color in [
                (f"riv25_{cohort}_baseline", "task_RIV25", "RIV 25 µg", "black"),
                (
                    f"riv25_{cohort}",
                    self.perpetrator_tasks[cohort],
                    f"RIV 25 µg + {self.perpetrator_labels[cohort]}",
                    self.colors[cohort],
                ),
            ]:
                plots[kp].add_data(
                    task=task,
                    xid="time",
                    yid="[Cve_riv]",
                    label=label,
                    color=color,
                )
                plots[kp].add_data(
                    dataset=key,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=color,
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Lenard2025, output_dir=Lenard2025.__name__)
