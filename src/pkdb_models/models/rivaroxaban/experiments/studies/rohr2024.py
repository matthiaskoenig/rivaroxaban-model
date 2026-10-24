"""Simulation experiment of the study Rohr2024."""

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


class Rohr2024(RivaroxabanSimulationExperiment):
    """Simulation experiment of Rohr2024.

    Microdose cocktail of apixaban 25 µg, rivaroxaban 25 µg and edoxaban 50 µg
    as oral solutions in eight healthy fasted subjects, at baseline and on the
    fifth day of ritonavir 100 mg twice daily ("RITx5"); plasma (Fig. 2,
    geometric means with the 95 % confidence interval) and the urine of 24 hr
    (Tab. 1, n = 7).

    Ritonavir inhibits CYP3A4 profoundly at 100 mg: the AUC of midazolam rose
    11- to 15-fold on days 2-5 in the same subjects, a residual CYP3A4 activity
    of about 5 % (`LI__f_cyp3a4`). It reduced the renal clearance of
    rivaroxaban from 49.4 to 34.8 ml/min (ratio 0.70), which the authors
    attribute to the inhibition of the renal efflux by P-gp and BCRP (and
    OAT3); the secretion is (RIVEX_f - 1) / RIVEX_f, about 87 % of the renal
    clearance, so the ratio 0.70 is a residual secretion of 65 %
    (`KI__f_pgp`). Observed: AUC ratio 1.87, Cmax ratio 1.27 (Tab. 1). The
    model has no intestinal P-gp or CYP3A4; the authors attribute the
    interaction to the systemic clearance (complete bioavailability of the
    microdose).
    """

    colors: ClassVar[dict[str, str]] = {
        "API25, RIV25, EDO50": "black",
        "API25, RIV25, EDO50, RITx5": "tab:red",
    }
    intervention_labels: ClassVar[dict[str, str]] = {
        "API25, RIV25, EDO50": "RIV25",
        "API25, RIV25, EDO50, RITx5": "RIV25 + ritonavir",
    }

    interventions: ClassVar[list[str]] = list(colors.keys())

    dose = 0.025  # [mg]

    #: the CYP3A4 and P-gp activities of the arms, see the class docstring
    inhibition: ClassVar[dict[str, dict[str, float]]] = {
        "API25, RIV25, EDO50": {"LI__f_cyp3a4": 1.0, "KI__f_pgp": 1.0},
        "API25, RIV25, EDO50, RITx5": {"LI__f_cyp3a4": 0.05, "KI__f_pgp": 0.65},
    }

    infos_pk: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Tab1A"]:
            df = self.load_dataframe(fig_id)
            # rivaroxaban only, not apixaban, edoxaban and the recovery
            df = df[
                (df.substance == "rivaroxaban")
                & df.measurement_type.isin(["concentration", "cumulative amount"])
            ]
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if "urine" not in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for intervention in self.interventions:
            simulations[self.to_sid(intervention)] = Simulation(
                time_unit="min",
                end=30 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "POSOLUTION_riv": Q(self.dose, "mg"),
                    **{
                        pid: Q(value, "dimensionless")
                        for pid, value in self.inhibition[intervention].items()
                    },
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            for sid, name in self.infos_pk.items():
                mappings[self.to_sid(f"fm_{name}_{intervention}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{intervention}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{intervention}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.RITONAVIR
                        if "RIT" in intervention
                        else Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {**self.figure_pk()}

    def figure_pk(self) -> dict[str, Figure]:
        """Figure of rivaroxaban in plasma and urine."""
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)

        for intervention in self.interventions:
            label = self.intervention_labels[intervention]
            for ks, (sid, name) in enumerate(self.infos_pk.items()):
                # simulation
                plots[ks].add_data(
                    task=self.to_sid(f"task_{intervention}"),
                    xid="time",
                    yid=sid,
                    label=label,
                    color=self.colors[intervention],
                )

                # data
                plots[ks].add_data(
                    dataset=self.to_sid(f"{name}_{intervention}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=label,
                    color=self.colors[intervention],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Rohr2024, output_dir=Rohr2024.__name__)
