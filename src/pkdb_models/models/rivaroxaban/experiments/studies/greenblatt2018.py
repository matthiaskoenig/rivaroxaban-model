from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)
from pkdb_models.models.rivaroxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, RivaroxabanMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rivaroxaban.helpers import run_experiments


class Greenblatt2018(RivaroxabanSimulationExperiment):
    """Simulation experiment of Greenblatt2018."""

    groups = ["NRF", "MRI"]
    colors = {
        "NRF": "black",
        "MRI": "#2ca25f",
    }
    bodyweights = {
        "NRF": 76.5,
        "MRI": 68.7,
    }
    group_labels = {
        "NRF": "Normal renal function",
        "MRI": "Mild renal impairment",
    }
    renal_functions = {
        "NRF": 105/105, # Normalized on NRI
        "MRI": 71/105,  # Normalized on NRI
    }
    infos = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}

        for fig_id in ["Fig2", "Fig5"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for group in self.groups:
            tcsims[f"po_{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=80 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[group], "kg"),
                        "KI__f_renal_function": Q_(self.renal_functions[group], "dimensionless"),
                        "PODOSE_riv": Q_(20, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            for k, sid in enumerate(self.infos):
                name = self.infos[sid]
                mappings[f"fm_po_{group}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_po_{group}", xid="time", yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY if group == "NRF" else Health.RENAL_IMPAIRMENT,
                        fasting=Fasting.FED,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_rows=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)

        for group in self.groups:
            for k, sid in enumerate(self.infos):
                name = self.infos[sid]

                # simulation
                plots[k].add_data(
                    task=f"task_po_{group}",
                    xid="time",
                    yid=sid,
                    label=self.group_labels[group],
                    color=self.colors[group],
                )

                # data
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.group_labels[group],
                    color=self.colors[group],
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Greenblatt2018.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Greenblatt2018, output_dir=Greenblatt2018.__name__)
