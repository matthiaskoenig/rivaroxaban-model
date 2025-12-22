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


class GenisNajera2022(RivaroxabanSimulationExperiment):
    """Simulation experiment of GenisNajera2022."""
    # reference and test formulation in fed conditions

    formulations = ["T20", "R20"]
    colors = {
        "R20": "black",
        "T20": "tab:blue",
    }
    bodyweight = 61.7

    formulation_labels = {
        "T20": "Test 20 mg",
        "R20": "Reference 20 mg",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_id = "Fig1"
        df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
        for label, df_label in df.groupby("label"):
            dset = DataSet.from_df(df_label, self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.riv)
            dsets[label] = dset

        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for formulation in self.formulations:
            tcsims[f"po_{formulation}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "PODOSE_riv": Q_(20, "mg"),
                        "GU__F_riv_abs": Q_(self.fasting_map["fed"], "dimensionless"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for formulation in self.formulations:
            mappings[f"fm_po_{formulation}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"rivaroxaban_{formulation}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self, task=f"task_po_{formulation}", xid="time", yid=f"[Cve_riv]",
                ),
                metadata=RivaroxabanMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.FED,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for formulation in self.formulations:

            # simulation
            plots[0].add_data(
                task=f"task_po_{formulation}",
                xid="time",
                yid=f"[Cve_riv]",
                label=f"{self.formulation_labels[formulation]}",
                color=self.colors[f"{formulation}"],
            )

            # data
            plots[0].add_data(
                dataset=f"rivaroxaban_{formulation}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{self.formulation_labels[formulation]}",
                color=self.colors[f"{formulation}"],
            )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / GenisNajera2022.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(GenisNajera2022, output_dir=GenisNajera2022.__name__)
