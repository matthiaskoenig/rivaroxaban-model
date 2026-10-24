"""Simulation experiment of the study Graff2007."""

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


class Graff2007(RivaroxabanSimulationExperiment):
    """Simulation experiment of Graff2007.

    12 healthy men, 8 of them took 5 mg and 30 mg rivaroxaban in a crossover
    after a fast of 10 hr, 4 took placebo; the form is not reported. Fig3 is
    the FXa inhibition and the PT relative to baseline (mean +- SD), Fig5 the
    geometric mean of the plasma concentration without an error, mapped as a
    mean. The placebo arm is shown, but not mapped: the model is exactly at
    its baseline without a dose, so the arm carries no information.
    """

    bodyweight = (66 + 90) / 2  # [kg] average bodyweight from range
    colors: ClassVar[dict[str, str]] = {
        "PLACEBO": "black",
        "RIV_5": "#F4A259",
        "RIV_30": "#D25100",
    }
    formulations: ClassVar[list[str]] = list(colors.keys())

    formulation_labels: ClassVar[dict[str, str]] = {
        "PLACEBO": "Placebo",
        "RIV_5": "Rivaroxaban 5 mg",
        "RIV_30": "Rivaroxaban 30 mg",
    }

    info_fig5: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
    }

    info_fig3: ClassVar[dict[str, str]] = {
        "Xa_inhibition": "factor Xa inhibition",
        "PT_ratio": "prothrombin time (change relative)",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3A", "Fig3B", "Fig5"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for formulation in self.formulations:
            if "5" in formulation:
                dose = 5
            elif "30" in formulation:
                dose = 30
            else:
                dose = 0

            tcsims[self.to_sid(f"{formulation}")] = Simulation(
                time_unit="min",
                end=40 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "PODOSE_riv": Q(dose, "mg"),
                },
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for formulation in self.formulations:
            if formulation == "PLACEBO":
                continue
            # pharmacokinetics
            for sid, name in self.info_fig5.items():
                mappings[self.to_sid(f"fm_{name}_{formulation}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{formulation}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{formulation}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.NR,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )

            # pharmacodynamics
            for sid, name in self.info_fig3.items():
                mappings[self.to_sid(f"fm_{name}_{formulation}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{formulation}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{formulation}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.NR,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_pk(),
            **self._figure_pd(),
        }

    def _figure_pk(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self, sid="Fig5", num_cols=2, name=self.__class__.__name__
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=1e-3
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False

        for formulation in self.formulations:
            if formulation == "PLACEBO":
                continue
            for sid, name in self.info_fig5.items():
                for plot in plots:
                    plot.add_data(
                        task=self.to_sid(f"task_{formulation}"),
                        xid="time",
                        yid=sid,
                        label=self.formulation_labels[formulation],
                        color=self.colors[formulation],
                    )
                    plot.add_data(
                        dataset=self.to_sid(f"{name}_{formulation}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=self.formulation_labels[formulation],
                        color=self.colors[formulation],
                    )

        return self.figures_by_sid(fig)

    def _figure_pd(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self, sid="Fig3", num_rows=2, name=self.__class__.__name__
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(
            self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"]
        )
        plots[1].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        for formulation in self.formulations:
            for k, (sid, name) in enumerate(self.info_fig3.items()):
                plots[k].add_data(
                    task=self.to_sid(f"task_{formulation}"),
                    xid="time",
                    yid=sid,
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{formulation}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Graff2007, output_dir=Graff2007.__name__)
