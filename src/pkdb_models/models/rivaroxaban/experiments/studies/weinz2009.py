"""Simulation experiment of the study Weinz2009."""

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
from sbmlsim.plot import Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)


class Weinz2009(RivaroxabanSimulationExperiment):
    """Simulation experiment of Weinz2009.

    Mass balance of a single oral solution of about 10 mg [14C]rivaroxaban in
    four healthy fasted men: unchanged rivaroxaban and the total radioactivity
    (rivaroxaban and its metabolites in rivaroxaban equivalents) in plasma
    (Fig. 2, geometric means with SD), the cumulative excretion of the
    radioactivity in urine and feces (Fig. 3, means with SD) and the unchanged
    rivaroxaban in the urine of 0-2 days (36.2 % of the dose) and in the feces
    of 0-4 days (7.3 %, Tab. 3).

    The SD of the plasma concentrations are read from error bars which the
    symbols of the log axis hide (SD 0 at 1 hr): the plasma is fitted without
    them.
    """

    dose = 10.0  # [mg]

    infos: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
        "[Cve_riv_rx]": "riv_rx",
        "Aurine_riv": "rivaroxaban_urine",
        "Afeces_riv": "rivaroxaban_feces",
        "Aurine_riv_rx": "riv_rx_urine",
        "Afeces_riv_rx": "riv_rx_feces",
        "Afeces_urine_riv_rx": "riv_rx_feces_urine",
    }
    tissues: ClassVar[dict[str, Tissue]] = {
        "rivaroxaban": Tissue.PLASMA,
        "riv_rx": Tissue.PLASMA,
        "rivaroxaban_urine": Tissue.URINE,
        "rivaroxaban_feces": Tissue.FECES,
        "riv_rx_urine": Tissue.URINE,
        "riv_rx_feces": Tissue.FECES,
        "riv_rx_feces_urine": Tissue.FECES_URINE,
    }
    labels_data: ClassVar[dict[str, str]] = {
        "rivaroxaban": "rivaroxaban",
        "riv_rx": "total radioactivity",
        "rivaroxaban_urine": "rivaroxaban urine",
        "rivaroxaban_feces": "rivaroxaban feces",
        "riv_rx_urine": "total urine",
        "riv_rx_feces": "total feces",
        "riv_rx_feces_urine": "total urine + feces",
    }
    colors: ClassVar[dict[str, str]] = {
        "rivaroxaban": "black",
        "riv_rx": "tab:blue",
        "rivaroxaban_urine": "tab:purple",
        "rivaroxaban_feces": "tab:brown",
        "riv_rx_urine": "tab:red",
        "riv_rx_feces": "tab:orange",
        "riv_rx_feces_urine": "tab:green",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            # the amounts, the recovery in percent of Fig. 3 is the same data
            df = df[df.measurement_type != "recovery"]
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset

        # unchanged rivaroxaban in urine and feces in percent of the dose
        df = self.load_dataframe("Tab3")
        df = df[df.substance == "rivaroxaban"]
        for tissue, df_tissue in df.groupby("tissue"):
            dset = DataSet.from_df(df_tissue, self.ureg)
            dset.unit_conversion("mean", Q(self.dose, "mg") / self.Mr.riv)
            dsets[self.to_sid(f"rivaroxaban_{tissue}")] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations["riv"] = Simulation(
            time_unit="min",
            end=8 * 24 * 60,
            preinit_changes={
                **self.default_changes(),
                "POSOLUTION_riv": Q(self.dose, "mg"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.infos.items():
            tissue = self.tissues[name]
            mappings[self.to_sid(f"fm_{name}")] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=self.to_sid(name),
                    xid="time",
                    yid="mean",
                    yid_sd=None if tissue == Tissue.PLASMA else "mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task="task_riv",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=tissue,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FASTED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        # Figure 2 (plasma, linear and log axis) and Figure 3 (excretion)
        fig = Figure(
            experiment=self,
            sid="Fig2_3",
            num_rows=1,
            num_cols=3,
            name=self.__class__.__name__,
            settings={"legend_scope": "row"},
        )
        plots = fig.create_plots(legend=True)
        for kp in [0, 1]:
            plots[kp].set_xaxis(self.label_time, unit=self.unit_time, min=-1, max=30)
            plots[kp].set_yaxis(
                self.label_riv_plasma,
                unit=self.unit_riv,
                scale="linear" if kp == 0 else "log",
                min=None if kp == 0 else 1e-3,
            )
        plots[2].set_xaxis(self.label_time, unit="day")
        plots[2].set_yaxis("excretion (14C)", unit=self.unit_riv_rx_feces)

        for sid, name in self.infos.items():
            kps = [0, 1] if self.tissues[name] == Tissue.PLASMA else [2]
            for kp in kps:
                plots[kp].add_data(
                    task="task_riv",
                    xid="time",
                    yid=sid,
                    label=self.labels_data[name],
                    color=self.colors[name],
                )
                plots[kp].add_data(
                    dataset=self.to_sid(name),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.labels_data[name],
                    color=self.colors[name],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Weinz2009, output_dir=Weinz2009.__name__)
