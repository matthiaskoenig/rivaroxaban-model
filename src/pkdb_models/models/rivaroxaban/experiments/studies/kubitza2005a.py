"""Simulation experiment of the study Kubitza2005a."""

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


class Kubitza2005a(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kubitza2005a.

    Single doses in healthy men after an overnight fast: tablets of 1.25 to
    80 mg (T) and the oral solution of 5 and 10 mg (S). Plasma of the tablets
    of 1.25, 5, 10, 20, 40 and 80 mg (Fig5A) and of the solutions with the 5
    and 10 mg tablets (Fig5B) as geometric means, factor Xa inhibition,
    prothrombin time and aPTT relative to the baseline (Fig3A-C) as medians;
    the plasma of the 15, 30 and 60 mg tablets is not shown in the
    publication.
    """

    formulations: ClassVar[list[str]] = [
        "T1.25",
        "T5",
        "T10",
        "T20",
        "T40",
        "T80",
        "S5",
        "S10",
    ]
    pd_formulations: ClassVar[list[str]] = [
        "P",
        "T1.25",
        "T5",
        "T10",
        "T20",
        "T40",
        "T80",
    ]
    bodyweight = 81.2

    colors: ClassVar[dict[str, str]] = {
        "T1.25": "#E6C229",
        "T5": "#FFD580",  # Light orange (peachy)
        "T10": "#FFB347",  # Soft orange
        "T20": "#FF7043",  # Reddish-orange
        "T40": "#A63C00",  # Darker brownish orange
        "T80": "#4C1F00",  # Very dark brown
        "S5": "tab:blue",
        "S10": "tab:purple",
        "P": "black",  # Placebo
    }

    formulation_labels: ClassVar[dict[str, str]] = {
        "P": "Placebo",
        "T1.25": "Rivaroxaban 1.25 mg",
        "T5": "Rivaroxaban 5 mg",
        "T10": "Rivaroxaban 10 mg",
        "T20": "Rivaroxaban 20 mg",
        "T40": "Rivaroxaban 40 mg",
        "T80": "Rivaroxaban 80 mg",
        "S5": "Rivaroxaban 5 mg solution",
        "S10": "Rivaroxaban 10 mg solution",
    }

    info_fig5: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
    }

    info_fig3: ClassVar[dict[str, str]] = {
        "Xa_inhibition": "factor Xa inhibition",
        "PT_ratio": "prothrombin time (change relative)",
        "aPTT_ratio": "aPTT (change relative)",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig5", "Fig3A", "Fig3B", "Fig3C"]

        for fig_id in fig_ids:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if df_label["measurement_type"].iloc[0] == "concentration":
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(f"{fig_id}_{label}")] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        all_formulations = list(dict.fromkeys(self.formulations + self.pd_formulations))

        for f in all_formulations:
            dose = 0.0 if f == "P" else float(f[1:])
            # the oral solution is the dissolved dose
            dose_id = "POSOLUTION_riv" if f.startswith("S") else "PODOSE_riv"
            tcsims[self.to_sid(f)] = Simulation(
                time_unit="min",
                end=50 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    dose_id: Q(dose, "mg"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                },
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        for f in self.formulations:
            for sid, name in self.info_fig5.items():
                dataset_id = f"Fig5_{name}_{f}"
                mappings[self.to_sid(f"fm_{dataset_id}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(dataset_id),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{f}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET
                        if f.startswith("T")
                        else ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )

        for f in self.pd_formulations:
            for fig_id, (sid, name) in zip(
                ["Fig3A", "Fig3B", "Fig3C"], self.info_fig3.items(), strict=True
            ):
                dataset_id = f"{fig_id}_{name}_{f}"
                mappings[self.to_sid(f"fm_{dataset_id}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(dataset_id),
                        xid="time",
                        yid="median",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{f}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        tablet_formulations = [f for f in self.formulations if f.startswith("T")]

        fig_pk1 = Figure(
            experiment=self,
            sid="Fig5A",
            name=self.__class__.__name__,
            settings={"legend_position": "outside"},
        )
        plots = fig_pk1.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=30), legend=True
        )

        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        for f in tablet_formulations:
            plots[0].add_data(
                task=self.to_sid(f"task_{f}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )
            plots[0].add_data(
                dataset=self.to_sid(f"Fig5_rivaroxaban_{f}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )

        fig_pd = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=3,
            name=self.__class__.__name__,
            settings={"legend_position": "outside"},
        )
        plots = fig_pd.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=30), legend=False
        )
        plots[0].legend = True

        for i, (sid, label) in enumerate(self.info_fig3.items()):
            plots[i].set_yaxis(self.labels[sid], unit=self.units[sid])
            for f in self.pd_formulations:
                plots[i].add_data(
                    task=self.to_sid(f"task_{f}"),
                    xid="time",
                    yid=sid,
                    label=self.formulation_labels[f],
                    color=self.colors[f],
                )
                plots[i].add_data(
                    dataset=self.to_sid(f"Fig3{chr(65 + i)}_{label}_{f}"),
                    xid="time",
                    yid="median",
                    count="count",
                    label=self.formulation_labels[f],
                    color=self.colors[f],
                )

        # Tablet and Solution
        fig_pk2 = Figure(
            experiment=self, sid="Fig5B", num_rows=2, name=self.__class__.__name__
        )
        plots = fig_pk2.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=30), legend=True
        )
        for kp in range(2):
            plots[kp].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[0].set_title("Oral solution")
        for f in ["S5", "S10"]:
            plots[0].add_data(
                task=self.to_sid(f"task_{f}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )
            plots[0].add_data(
                dataset=self.to_sid(f"Fig5_rivaroxaban_{f}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )

        plots[1].set_title("Tablet")
        for f in ["T5", "T10"]:
            plots[1].add_data(
                task=self.to_sid(f"task_{f}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )
            plots[1].add_data(
                dataset=self.to_sid(f"Fig5_rivaroxaban_{f}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.formulation_labels[f],
                color=self.colors[f],
            )

        return self.figures_by_sid(fig_pk1, fig_pd, fig_pk2)


if __name__ == "__main__":
    MODEL.run_experiments(Kubitza2005a, output_dir=Kubitza2005a.__name__)
