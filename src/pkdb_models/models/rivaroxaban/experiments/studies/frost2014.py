from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)
from pkdb_models.models.rivaroxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, RivaroxabanMappingMetaData

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.rivaroxaban.helpers import run_experiments


class Frost2014(RivaroxabanSimulationExperiment):
    """Simulation experiment of Frost2014.

    Rivaroxaban once daily, 10 mg
    """

    bodyweight = 75.9 # [kg]

    colors = {
        "mRIV10": "black"
    }
    interventions = list(colors.keys())
    doses = {
        "mRIV10": 10
    }
    infos_pk = {
        "[Cve_riv]": "rivaroxaban"
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset

        # console.print(dsets)

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            dose = self.doses[intervention]
            tc0 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=100,
                changes={
                    **self.default_changes(),
                    "BW": Q_(self.bodyweight, "kg"),
                    "PODOSE_riv": Q_(dose, "mg"),
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=100,
                changes={
                    "PODOSE_riv": Q_(dose, "mg"),
                },
            )
            tc2 = Timecourse(
                start=0,
                end=75 * 60,  # [min]
                steps=1000,
                changes={
                    "PODOSE_riv": Q_(dose, "mg"),
                },
            )

            tcsims[intervention] = TimecourseSim(
                [tc0] + [tc1 for _ in range(2)] + [tc2],
                time_offset=-3 * 24 * 60
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            # pharmacokinetics
            for k, sid in enumerate(self.infos_pk.keys()):

                name = self.infos_pk[sid]
                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing= Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting= Fasting.NR,
                    )
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig2",
            name=f"{self.__class__.__name__} (healthy)"
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for intervention in self.interventions:
            for k, sid in enumerate(self.infos_pk.keys()):
                name = self.infos_pk[sid]
                #simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=intervention,
                    color=self.colors[intervention],
                )
                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd=None,
                    count="count",
                    label=intervention,
                    color=self.colors[intervention],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Frost2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Frost2014, output_dir=Frost2014.__name__)
