"""Simulation experiment of the study Kubitza2006."""

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


class Kubitza2006(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kubitza2006.

    Food and formulation: 10 mg (2 x 5 mg tablets) fasted or after a high-fat
    breakfast (Fig2 plasma, Fig3 prothrombin time, Fig4 factor Xa inhibition)
    and 20 mg as 4 x 5 mg fasted, 1 x 20 mg fasted or fed (Fig5 plasma) in
    healthy men. The ranitidine and antacid arms have no time courses.
    """

    groups: ClassVar[list[str]] = [
        "RIV10_fed",
        "RIV10_fasted",
        "RIV20_fed",
        "RIV20_fasted",
        "RIV4x5_fasted",
    ]

    group_labels: ClassVar[dict[str, str]] = {
        "RIV10_fed": "Rivaroxaban 10 mg (fed)",
        "RIV10_fasted": "Rivaroxaban 10 mg (fasted)",
        "RIV20_fed": "Rivaroxaban 20 mg (fed)",
        "RIV20_fasted": "Rivaroxaban 20 mg (fasted)",
        "RIV4x5_fasted": "Rivaroxaban 4 x 5 mg (fasted)",
    }

    bodyweights: ClassVar[dict[str, float]] = {
        "RIV10_fed": 85,
        "RIV10_fasted": 85,
        "RIV20_fed": 86,
        "RIV20_fasted": 86,
        "RIV4x5_fasted": 86,
    }

    doses: ClassVar[dict[str, float]] = {
        "RIV10_fed": 10,
        "RIV10_fasted": 10,
        "RIV20_fed": 20,
        "RIV20_fasted": 20,
        "RIV4x5_fasted": 20,
    }

    colors: ClassVar[dict[str, str]] = {
        "RIV10_fed": "#A35AC7",
        "RIV10_fasted": "#DDB8EA",
        "RIV20_fed": "#A35AC7",
        "RIV20_fasted": "#DDB8EA",
        "RIV4x5_fasted": "tab:gray",
    }

    info_fig2_5: ClassVar[dict[str, str]] = {"[Cve_riv]": "rivaroxaban"}

    info_fig3_4: ClassVar[dict[str, str]] = {
        "PT_ratio": "prothrombin time (change relative)",
        "Xa_inhibition": "factor Xa inhibition",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig2", "Fig3", "Fig4", "Fig5"]

        for fig_id in fig_ids:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if df_label["measurement_type"].iloc[0] == "concentration":
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        for group in self.groups:
            dose = self.doses[group]
            bw = self.bodyweights[group]
            fed = self.fasting_map["fed" if "fed" in group else "fasted"]

            tcsims[self.to_sid(group)] = Simulation(
                time_unit="min",
                end=30 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(bw, "kg"),
                    "PODOSE_riv": Q(dose, "mg"),
                    "GU__fed": Q(fed, "dimensionless"),
                },
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        # PK mapping for all groups
        for group in self.groups:
            fed = "fed" in group
            for sid, name in self.info_fig2_5.items():
                mappings[self.to_sid(f"fm_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{group}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if fed else Fasting.FASTED,
                    ),
                )

        # PD mappings only for the correct groups
        pd_groups = ["RIV10_fed", "RIV10_fasted"]

        for group in pd_groups:
            fed = "fed" in group
            for sid, name in self.info_fig3_4.items():
                mappings[self.to_sid(f"fm_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{group}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if fed else Fasting.FASTED,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_concentration(),
            **self._figure_pd(),
        }

    def _figure_concentration(self) -> dict[str, Figure]:
        fig1 = Figure(
            experiment=self, sid="Fig_conc_10mg", name=self.__class__.__name__
        )

        fig5 = Figure(
            experiment=self, sid="Fig_conc_20mg", name=self.__class__.__name__
        )

        plot_fig1 = fig1.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plot_fig1[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        plot_fig5 = fig5.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plot_fig5[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for group in ["RIV10_fed", "RIV10_fasted"]:
            plot_fig1[0].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plot_fig1[0].add_data(
                dataset=self.to_sid(f"rivaroxaban_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

        for group in ["RIV20_fed", "RIV20_fasted", "RIV4x5_fasted"]:
            plot_fig5[0].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plot_fig5[0].add_data(
                dataset=self.to_sid(f"rivaroxaban_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
                marker="^" if "4x5" in group else "s",
            )

        return self.figures_by_sid(fig1, fig5)

    def _figure_pd(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_PD",
            num_cols=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        plots[0].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        plots[1].set_yaxis(
            self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"]
        )

        for group in ["RIV10_fed", "RIV10_fasted"]:
            plots[0].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="PT_ratio",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plots[0].add_data(
                dataset=self.to_sid(f"prothrombin time (change relative)_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

            plots[1].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="Xa_inhibition",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plots[1].add_data(
                dataset=self.to_sid(f"factor Xa inhibition_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kubitza2006, output_dir=Kubitza2006.__name__)
