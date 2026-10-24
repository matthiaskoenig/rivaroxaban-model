"""Simulation experiment of the study Greenblatt2018."""

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


class Greenblatt2018(RivaroxabanSimulationExperiment):
    """Simulation experiment of Greenblatt2018.

    Subjects with a normal renal function (NRF, creatinine clearance 105
    ml/min, n = 13) and a mild renal impairment (MRI, 71 ml/min, n = 14), mean
    age 59 yr, a single 20 mg tablet of Xarelto right after a standard
    breakfast (fed) on day 1 and on day 15. Verapamil extended release once
    daily on days 8 to 17 (120 mg, 240 mg, then 360 mg), so the dose of day 15
    is taken with verapamil at steady state (n = 10 and 11). Fig2 is the
    plasma (mean +- SD), Fig5 the cumulative urinary excretion (mean +- SD).
    """

    groups: ClassVar[list[str]] = ["NRF", "MRI"]
    #: the day 1 dose alone and the day 15 dose with verapamil
    interventions: ClassVar[list[str]] = ["", "VER"]
    colors: ClassVar[dict[str, str]] = {
        "NRF": "black",
        "MRI": "#2ca25f",
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "NRF": 76.5,  # kg
        "MRI": 68.7,  # kg
    }
    group_labels: ClassVar[dict[str, str]] = {
        "NRF": "normal renal function",
        "MRI": "mild renal impairment",
    }
    intervention_labels: ClassVar[dict[str, str]] = {
        "": "",
        "VER": " + verapamil",
    }
    linestyles: ClassVar[dict[str, str]] = {
        "": "-",
        "VER": "--",
    }
    markers: ClassVar[dict[str, str]] = {
        "": "s",
        "VER": "^",
    }
    #: the renal function relative to the GFR of 101 ml/min, the estimated
    #: creatinine clearance of the groups (Tab1)
    renal_functions: ClassVar[dict[str, float]] = {
        "NRF": 105 / 101,
        "MRI": 71 / 101,
    }
    #: verapamil 360 mg (extended release) at steady state, not fitted.
    #: Verapamil and norverapamil are mechanism-based inhibitors of CYP3A4,
    #: the AUC of oral midazolam rises 2.9-3.8x (Backman1994), i.e. an
    #: activity of 0.2-0.3 for a midazolam fraction metabolized by CYP3A4 of
    #: 0.93; 0.3 is taken. The renal secretion is not inhibited: verapamil is a
    #: weak P-gp inhibitor at its unbound plasma concentrations. Observed (Tab2,
    #: Fig5): AUC ratio 1.43 (NRF) and 1.39 (MRI), urinary excretion 1.55x and
    #: 1.37x; the model gives about 1.2x for both, the remainder is an
    #: increase of the absorption by the intestinal P-gp/BCRP, which the model
    #: does not have.
    verapamil_changes: ClassVar[dict[str, float]] = {
        "LI__f_cyp3a4": 0.3,
        "KI__f_pgp": 1.0,
    }
    infos: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
    }

    @staticmethod
    def _arm(group: str, intervention: str) -> str:
        """Key of an arm, the group and the intervention."""
        return f"{group}_{intervention}" if intervention else group

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig5"]:
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
        for group in self.groups:
            for intervention in self.interventions:
                changes = {
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[group], "kg"),
                    "KI__f_renal_function": Q(
                        self.renal_functions[group], "dimensionless"
                    ),
                    "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
                    "PODOSE_riv": Q(20, "mg"),
                }
                if intervention == "VER":
                    for sid, value in self.verapamil_changes.items():
                        changes[sid] = Q(value, "dimensionless")
                arm = self._arm(group, intervention)
                tcsims[self.to_sid(f"po_{arm}")] = Simulation(
                    time_unit="min",
                    end=80 * 60,
                    preinit_changes=changes,
                )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            for intervention in self.interventions:
                arm = self._arm(group, intervention)
                for sid, name in self.infos.items():
                    mappings[self.to_sid(f"fm_po_{arm}_{name}")] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=self.to_sid(f"{name}_{arm}"),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=self.to_sid(f"task_po_{arm}"),
                            xid="time",
                            yid=sid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY
                            if group == "NRF"
                            else Health.RENAL_IMPAIRMENT,
                            fasting=Fasting.FED,
                            coadministration=Coadministration.VERAPAMIL
                            if intervention == "VER"
                            else Coadministration.NONE,
                        ),
                    )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_5",
            num_rows=2,
            num_cols=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=1e-3
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False
        plots[2].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
        plots[3].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
        # the panels of the observables: the plasma linear and log, the urine
        # of day 1 and of day 15 with verapamil
        panels: dict[tuple[str, str], list[int]] = {
            ("[Cve_riv]", ""): [0, 1],
            ("[Cve_riv]", "VER"): [0, 1],
            ("Aurine_riv", ""): [2],
            ("Aurine_riv", "VER"): [3],
        }

        for group in self.groups:
            for intervention in self.interventions:
                arm = self._arm(group, intervention)
                label = (
                    f"{self.group_labels[group]}"
                    f"{self.intervention_labels[intervention]}"
                )
                for sid, name in self.infos.items():
                    for k in panels[(sid, intervention)]:
                        # simulation
                        plots[k].add_data(
                            task=self.to_sid(f"task_po_{arm}"),
                            xid="time",
                            yid=sid,
                            label=label,
                            color=self.colors[group],
                            linestyle=self.linestyles[intervention],
                        )

                        # data
                        plots[k].add_data(
                            dataset=self.to_sid(f"{name}_{arm}"),
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                            label=label,
                            color=self.colors[group],
                            marker=self.markers[intervention],
                        )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Greenblatt2018, output_dir=Greenblatt2018.__name__)
