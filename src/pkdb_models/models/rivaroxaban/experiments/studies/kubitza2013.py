"""Simulation experiment of the study Kubitza2013."""

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


class Kubitza2013(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kubitza2013.

    Single dose of rivaroxaban 10 mg (2 x 5 mg tablets) after an overnight fast
    in healthy subjects (n = 16) and in subjects with mild (Child-Pugh A) or
    moderate (Child-Pugh B) hepatic impairment (n = 8 each). Plasma (Fig1),
    factor Xa inhibition and prothrombin time relative to the baseline (Fig2)
    as mean and SD. The renal function is the creatinine clearance of Table 1
    relative to the GFR of 101 ml/min.
    """

    groups: ClassVar[list[str]] = ["healthy", "mildHI", "moderateHI"]
    colors: ClassVar[dict[str, str]] = {
        "healthy": "black",
        "mildHI": "#74A9CF",
        "moderateHI": "#2B8CBE",
    }
    group_labels: ClassVar[dict[str, str]] = {
        "healthy": "Healthy subjects",
        "mildHI": "Mild hepatic impairment",
        "moderateHI": "Moderate hepatic impairment",
    }
    bodyweights: ClassVar[dict[str, float]] = {
        "healthy": 76,
        "mildHI": 81.1,
        "moderateHI": 73.3,
    }
    #: creatinine clearance [ml/min] (Table 1)
    crcl: ClassVar[dict[str, float]] = {
        "healthy": 134.5,
        "mildHI": 101.3,
        "moderateHI": 130.8,
    }
    hepatic_functions: ClassVar[dict[str, str]] = {
        "healthy": "Control",
        "mildHI": "Mild cirrhosis",
        "moderateHI": "Moderate cirrhosis",
    }

    info_figs: ClassVar[dict[str, dict[str, str]]] = {
        "Fig1": {"[Cve_riv]": "rivaroxaban"},
        "Fig2A": {"Xa_inhibition": "factor Xa inhibition"},
        "Fig2B": {"PT_ratio": "prothrombin time (change relative)"},
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2A", "Fig2B"]:
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
        simulations = {}
        for group in self.groups:
            hf_value = self.cirrhosis_map[self.hepatic_functions[group]]
            simulations[self.to_sid(group)] = Simulation(
                time_unit="min",
                end=100 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[group], "kg"),
                    "PODOSE_riv": Q(10, "mg"),
                    "f_cirrhosis": Q(hf_value, "dimensionless"),
                    "KI__f_renal_function": Q(self.crcl[group] / 101, "dimensionless"),
                    "GU__fed": Q(self.fasting_map["fasted"], "dimensionless"),
                },
            )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for fig_id, variables in self.info_figs.items():
            for yid, label in variables.items():
                for group in self.groups:
                    dataset_id = f"{fig_id}_{label}_{group}"
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
                            task=self.to_sid(f"task_{group}"),
                            xid="time",
                            yid=yid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY
                            if group == "healthy"
                            else Health.HEPATIC_IMPAIRMENT,
                            fasting=Fasting.FASTED,
                        ),
                    )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        # PK figure: Fig1 rivaroxaban concentration
        fig_pk = Figure(experiment=self, sid="Fig1", name=self.__class__.__name__)
        plots_pk = fig_pk.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=50),
            legend=True,
        )
        plots_pk[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        for group in self.groups:
            plots_pk[0].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="[Cve_riv]",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plots_pk[0].add_data(
                dataset=self.to_sid(f"Fig1_rivaroxaban_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

        # PD figure: Fig2A and Fig2B side by side
        fig_pd = Figure(
            experiment=self, sid="Fig2", num_cols=2, name=self.__class__.__name__
        )
        plots_pd = fig_pd.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, max=75),
            legend=True,
        )

        # Fig2A: Xa inhibition (left plot)
        plots_pd[0].set_yaxis(
            self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"]
        )
        for group in self.groups:
            plots_pd[0].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="Xa_inhibition",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plots_pd[0].add_data(
                dataset=self.to_sid(f"Fig2A_factor Xa inhibition_{group}"),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

        # Fig2B: prothrombin time (right plot)
        plots_pd[1].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        for group in self.groups:
            plots_pd[1].add_data(
                task=self.to_sid(f"task_{group}"),
                xid="time",
                yid="PT_ratio",
                label=self.group_labels[group],
                color=self.colors[group],
            )
            plots_pd[1].add_data(
                dataset=self.to_sid(
                    f"Fig2B_prothrombin time (change relative)_{group}"
                ),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.group_labels[group],
                color=self.colors[group],
            )

        return self.figures_by_sid(fig_pk, fig_pd)


if __name__ == "__main__":
    MODEL.run_experiments(Kubitza2013, output_dir=Kubitza2013.__name__)
