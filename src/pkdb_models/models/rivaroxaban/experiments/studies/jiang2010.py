"""Simulation experiment of the study Jiang2010."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
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


class Jiang2010(RivaroxabanSimulationExperiment):
    """Simulation experiment of Jiang2010.

    Single oral doses of 5 to 40 mg in healthy elderly Chinese subjects, taken
    within 5 min after a standardized breakfast (fed). Plasma (Fig1, geometric
    means), urinary excretion 0-72 hr (Table 2), factor Xa inhibition (Fig2)
    and prothrombin time (Fig3A) as medians. The subsets of male and female
    subjects (Fig2B, Fig3B) are not simulated.
    """

    formulations: ClassVar[list[str]] = ["P", "T5", "T10", "T20", "T30", "T40"]
    colors: ClassVar[dict[str, str]] = {
        "P": "black",
        "T5": "#FFDAB9",
        "T10": "#FFB07C",
        "T20": "#FF874C",
        "T30": "#D25100",
        "T40": "#892E0E",
    }

    formulation_labels: ClassVar[dict[str, str]] = {
        "P": "Placebo",
        "T5": "Rivaroxaban 5 mg",
        "T10": "Rivaroxaban 10 mg",
        "T20": "Rivaroxaban 20 mg",
        "T30": "Rivaroxaban 30 mg",
        "T40": "Rivaroxaban 40 mg",
    }

    bodyweights: ClassVar[dict[str, float]] = {
        "P": 62.1,
        "T5": 62.5,
        "T10": 59.2,
        "T20": 60.5,
        "T30": 62.2,
        "T40": 63.5,
    }

    #: observable: (figure, label of the data, column of the central value)
    info_pk: ClassVar[dict[str, tuple[str, str, str]]] = {
        "[Cve_riv]": ("Fig1", "rivaroxaban", "mean"),
        "Aurine_riv": ("Tab2A", "rivaroxaban_urine", "mean"),
    }

    #: the pharmacodynamics are medians
    info_pd: ClassVar[dict[str, tuple[str, str, str]]] = {
        "Xa_inhibition": ("Fig2", "factor Xa inhibition", "median"),
        "PT_ratio": ("Fig3", "prothrombin time (change relative)", "median"),
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Fig3", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                measurement = df_label["measurement_type"].iloc[0]
                if measurement == "recovery":
                    # the same excretion in percent of the dose
                    continue
                dset = DataSet.from_df(df_label, self.ureg)
                if measurement in {"concentration", "cumulative amount"}:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(f"{fig_id}_{label}")] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for formulation in self.formulations:
            dose = 0 if formulation == "P" else float(formulation[1:])

            tcsims[self.to_sid(f"{formulation}")] = Simulation(
                time_unit="min",
                end=100 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[formulation], "kg"),
                    "PODOSE_riv": Q(dose, "mg"),
                    "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
                },
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for formulation in self.formulations:
            infos = {**self.info_pk, **self.info_pd}
            for sid, (fig_id, name, yid) in infos.items():
                if formulation == "P" and sid in self.info_pk:
                    continue
                urine = sid == "Aurine_riv"
                dataset = f"{fig_id}_{name}_{formulation}"
                mappings[self.to_sid(f"fm_{dataset}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(dataset),
                        xid="time",
                        yid=yid,
                        yid_sd="mean_sd" if urine else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{formulation}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if urine else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure(sid="Fig1_Tab2A", info=self.info_pk, num_cols=2),
            **self._figure(sid="Fig2_3", info=self.info_pd, num_cols=2),
        }

    def _figure(
        self, sid: str, info: dict[str, tuple[str, str, str]], num_cols: int
    ) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid=sid,
            num_cols=num_cols,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        for k, (yid, (fig_id, name, central)) in enumerate(info.items()):
            plots[k].set_yaxis(
                self.labels[yid].replace("\n", " "), unit=self.units[yid]
            )
            for formulation in self.formulations:
                plots[k].add_data(
                    task=self.to_sid(f"task_{formulation}"),
                    xid="time",
                    yid=yid,
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
                if formulation == "P" and yid in self.info_pk:
                    continue
                plots[k].add_data(
                    dataset=self.to_sid(f"{fig_id}_{name}_{formulation}"),
                    xid="time",
                    yid=central,
                    yid_sd="mean_sd" if yid == "Aurine_riv" else None,
                    count="count",
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                    linestyle="" if yid == "Aurine_riv" else "--",
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Jiang2010, output_dir=Jiang2010.__name__)
