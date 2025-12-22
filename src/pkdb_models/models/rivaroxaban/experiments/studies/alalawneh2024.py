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


class Alalawneh2024(RivaroxabanSimulationExperiment):
    """Simulation experiment of Alalawneh2024."""

    groups = ["non-obese", "obese", "morbidly-obese"]
    bodyweights = {
        "non-obese": 67.4,  # kg
        "obese": 111.5,  # kg
        "morbidly-obese": 135, # kg > 120 kg
    }
    colors = {
        "non-obese": "black",
        "obese": "tab:blue",
        "morbidly-obese": "tab:orange",
    }
    info_fig3 = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
        "KI__RIVEX": "rivaroxaban_excretion_urine",
    }
    info_fig4 = {
        "PT": "prothrombin time",
        "aPTT": "aPTT",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig2", "Fig3", "Fig4a", "Fig4b", "Tab4"]

        for fig_id in fig_ids:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset

        # print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for group in self.groups:
            tcsims[f"{group}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[group], "kg"),
                        "PODOSE_riv": Q_(20, "mg"),
                    },
                )]
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            # pharmacokinetics
            for k, sid in enumerate(self.info_fig3.keys()):
                if k != 0 and "morbid" in group:
                    continue
                name = self.info_fig3[sid]
                mappings[f"fm_{name}_{group}"] = FitMapping(
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
                        self,
                        task=f"task_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )

            # pharmacodynamics
            if "morbid" in group:
                continue
            for k, sid in enumerate(self.info_fig4.keys()):
                name = self.info_fig4[sid]
                mappings[f"fm_{name}_{group}"] = FitMapping(
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
                        self,
                        task=f"task_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )


        return mappings

    def figures(self) -> Dict[str, Figure]:

        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig2_3_Tab4",
            num_rows=3,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
        plots[2].set_yaxis(self.labels["KI__RIVEX"], unit=self.units["KI__RIVEX"])

        for group in self.groups:
            for k, sid in enumerate(self.info_fig3.keys()):
                name = self.info_fig3[sid]
                label = f"{name.replace('_', ' ')} {group}"
                if k != 0 and "morbid" in group:
                    continue

                plots[k].add_data(
                    task=f"task_{group}",
                    xid="time",
                    yid=sid,
                    label=group,
                    color=self.colors[group],
                )
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=group,
                    color=self.colors[group],
                )

        return {fig.sid: fig}


    def figure_pd(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4",
            num_rows=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.labels["PT"], unit=self.units["PT"])
        plots[1].set_yaxis(self.labels["aPTT"], unit=self.units["aPTT"])
        for group in self.groups:
            if "morbid" in group:
                continue
            for k, sid in enumerate(self.info_fig4.keys()):
                name = self.info_fig4[sid]
                label = f"{sid} {group}"
                plots[k].add_data(
                    task=f"task_{group}",
                    xid="time",
                    yid=sid,
                    label=group,
                    color=self.colors[group],
                )
                plots[k].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=group,
                    color=self.colors[group],
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Alalawneh2024.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Alalawneh2024, output_dir=Alalawneh2024.__name__)
