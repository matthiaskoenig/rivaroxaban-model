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


class Graff2007(RivaroxabanSimulationExperiment):
    """Simulation experiment of Graff2007."""


    bodyweight = (66 + 90)/2  # [kg] average bodyweight from range
    colors = {
        "PLACEBO": "black",
        "RIV_5": "#FFDAB9",
        "RIV_30": "#D25100",
    }
    formulations = list(colors.keys())

    formulation_labels = {
    "PLACEBO": "Placebo",
    "RIV_5": "Rivaroxaban 5 mg",
    "RIV_30": "Rivaroxaban 30 mg",
}


    info_fig5 = {
        "[Cve_riv]": "rivaroxaban",
    }

    info_fig3 = {
        "Xa_inhibition": "factor Xa inhibition",
        "PT_ratio": "prothrombin time (change relative)",
    }


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig3A", "Fig3B", "Fig5"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for formulation in self.formulations:
            if "5" in formulation:
                dose = 5
            elif "30" in formulation:
                dose = 30
            else:
                dose = 0

            tcsims[f"{formulation}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=40 * 60,  # [min]
                    steps=1000,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "PODOSE_riv": Q_(dose, "mg"),
                    },
                )]
            )

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for formulation in self.formulations:
            # pharmacokinetics
            if formulation != "PLACEBO":
                for k, sid in enumerate(self.info_fig5.keys()):
                    name = self.info_fig5[sid]
                    mappings[f"fm_{name}_{formulation}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{formulation}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=f"task_{formulation}",
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

            # pharmacodynamics
            for k, sid in enumerate(self.info_fig3.keys()):
                name = self.info_fig3[sid]
                mappings[f"fm_{name}_{formulation}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{formulation}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{formulation}",
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
            sid="Fig5",
            num_rows=1,
            name=self.__class__.__name__
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for formulation in self.formulations:

            for k, sid in enumerate(self.info_fig5.keys()):
                name = self.info_fig5[sid]

                plots[k].add_data(
                    task=f"task_{formulation}",
                    xid="time",
                    yid=sid,
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
                if formulation != "PLACEBO":
                    plots[k].add_data(
                        dataset=f"{name}_{formulation}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=self.formulation_labels[formulation],
                        color=self.colors[formulation],
                    )

        return {fig.sid: fig}


    def figure_pd(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=2,
            name=self.__class__.__name__
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"])
        plots[1].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        for formulation in self.formulations:
            for k, sid in enumerate(self.info_fig3.keys()):
                name = self.info_fig3[sid]
                label = f"{sid} {formulation}"
                plots[k].add_data(
                    task=f"task_{formulation}",
                    xid="time",
                    yid=sid,
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
                plots[k].add_data(
                    dataset=f"{name}_{formulation}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )
        return {fig.sid: fig}


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Graff2007.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Graff2007, output_dir=Graff2007.__name__)
