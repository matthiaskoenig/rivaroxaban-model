from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment
)
from pkdb_models.models.rivaroxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, Health, \
    Fasting, RivaroxabanMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rivaroxaban.helpers import run_experiments


class Lenard2025(RivaroxabanSimulationExperiment):
    """Simulation experiment of Lenard2025."""

    colors = {
        "API25, RIV25, EDO50_1": "black",
        # "API25, RIV25, EDO50, CAR": "tab:blue",
        "API25, RIV25, EDO50_2": "black",
        # "API25, RIV25, EDO50, GAB": "tab:orange",
        "API25, RIV25, EDO50_3": "black",
        # "API25, RIV25, EDO50, PRE": "tab:red"
    }
    labels = {
        "API25, RIV25, EDO50_1": "RIV25 (1)",
        # "API25, RIV25, EDO50, CAR": "RIV25, CAR",
        "API25, RIV25, EDO50_2": "RIV25 (2)",
        # "API25, RIV25, EDO50, GAB": "RIV25, GAB",
        "API25, RIV25, EDO50_3": "RIV25 (3)",
        # "API25, RIV25, EDO50, PRE": "RIV25, PRE",
    }
    interventions =  list(colors.keys())
    dose = 0.025
    infos_pk = {
        "[Cve_riv]": "rivaroxaban"
    }


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3", "Fig4"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("rivaroxaban"):
                        dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            tcsims[f"{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "PODOSE_riv": Q_(self.dose, "mg")
                    },
                )])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            # PK
            for ks, sid in enumerate(self.infos_pk):
                name = self.infos_pk[sid]
                coadministration = Coadministration.NONE
                if "CAR" in intervention:
                    coadministration = Coadministration.CARBAMAZEPINE
                elif "GAB" in intervention:
                    coadministration = Coadministration.GABAPENTIN
                elif "PRE" in intervention:
                    coadministration = Coadministration.PREGABALIN
                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_{intervention}", xid="time", yid=sid,
                    ),
                    metadata= RivaroxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=coadministration,

                    ),
                )

        return mappings


    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pk):
                    name = self.infos_pk[sid]
                    # simulation
                    plots[ks].add_data(
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                        label=self.labels[intervention],
                        color=self.colors[f"{intervention}"],
                    )
                    # data
                    plots[ks].add_data(
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=self.labels[intervention],
                        color=self.colors[intervention],
                    )

        return {
            fig.sid: fig
        }


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Lenard2025.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Lenard2025, output_dir=Lenard2025.__name__)