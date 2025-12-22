from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import RivaroxabanSimulationExperiment
from pkdb_models.models.rivaroxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, RivaroxabanMappingMetaData
from pkdb_models.models.rivaroxaban.helpers import run_experiments

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim


class Kubitza2013(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kubitza2013."""

    groups = ["healthy", "mildHI", "moderateHI"]
    colors = {
        "healthy": "black",
        "mildHI": "#74A9CF",
        "moderateHI": "#2B8CBE",
    }
    group_labels = {
        "healthy": "Healthy subjects",
        "mildHI": "Mild hepatic impairment",
        "moderateHI": "Moderate hepatic impairment",
    }
    bodyweights = {
        "healthy": 76,
        "mildHI": 81.1,
        "moderateHI": 73.3,
    }
    hepatic_functions = {
        "healthy": "Control",
        "mildHI": "Mild cirrhosis",
        "moderateHI": "Moderate cirrhosis",
    }

    info_figs = {
        "Fig1": {
            "[Cve_riv]": "rivaroxaban"
        },
        "Fig2A": {
            "Xa_inhibition": "factor Xa inhibition"
        },
        "Fig2B": {
            "PT_ratio": "prothrombin time (change relative)"
        }
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig1", "Fig2A", "Fig2B"]

        for fig_id in fig_ids:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[f"{fig_id}_{label}"] = dset

        # print(dsets.keys())
        return dsets


    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for group in self.groups:
            hf = self.hepatic_functions[group]
            hf_value = RivaroxabanSimulationExperiment.cirrhosis_map[hf]
            tcsims[group] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=100 * 60,  # min
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[group], "kg"),
                        "PODOSE_riv": Q_(10, "mg"),
                        "f_cirrhosis": Q_(hf_value, "dimensionless"),
                    }
                )
            ])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for fig_id, variables in self.info_figs.items():
            for yid, label in variables.items():
                for group in self.groups:
                    dataset_id = f"{fig_id}_{label}_{group}"
                    # print(dataset_id)
                    mappings[f"fm_{dataset_id}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=dataset_id,
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=f"task_{group}",
                            xid="time",
                            yid=yid,
                        ),
                        metadata=RivaroxabanMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY if group == "healthy" else Health.HEPATIC_IMPAIRMENT,
                            fasting=Fasting.FASTED,
                        ),
                    )
        return mappings


    def figures(self) -> Dict[str, Figure]:
        figures = {}

        # PK figure: Fig1 rivaroxaban concentration
        fig_pk = Figure(
            experiment=self,
            sid="Fig1",
            name=self.__class__.__name__
        )
        plots_pk = fig_pk.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots_pk[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        for group in self.groups:
            plots_pk[0].add_data(
                task=f"task_{group}",
                xid="time",
                yid="[Cve_riv]",
                label = self.group_labels[group],
                color=self.colors[group],
            )
            plots_pk[0].add_data(
                dataset=f"Fig1_rivaroxaban_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label = self.group_labels[group],
                color=self.colors[group],
            )
        figures[fig_pk.sid] = fig_pk

        # PD figure: Fig2A and Fig2B combined, 2 subplots
        fig_pd = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=2,
            name=self.__class__.__name__
        )
        plots_pd = fig_pd.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        # PD figure: Fig2A and Fig2B combined, now horizontally aligned
        fig_pd = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=self.__class__.__name__
        )
        plots_pd = fig_pd.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        # Fig2A: Xa inhibition (left plot)
        plots_pd[0].set_yaxis(self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"])
        for group in self.groups:
            plots_pd[0].add_data(
                task=f"task_{group}",
                xid="time",
                yid="Xa_inhibition",
                label = self.group_labels[group],
                color=self.colors[group],
            )
            plots_pd[0].add_data(
                dataset=f"Fig2A_factor Xa inhibition_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label = self.group_labels[group],
                color=self.colors[group],
            )

        # Fig2B: prothrombin time (right plot)
        plots_pd[1].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        for group in self.groups:
            plots_pd[1].add_data(
                task=f"task_{group}",
                xid="time",
                yid="PT_ratio",
                label = self.group_labels[group],
                color=self.colors[group],
            )
            plots_pd[1].add_data(
                dataset=f"Fig2B_prothrombin time (change relative)_{group}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label = self.group_labels[group],
                color=self.colors[group],
            )

        figures[fig_pd.sid] = fig_pd
        return figures


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Kubitza2013.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kubitza2013, output_dir=Kubitza2013.__name__)
