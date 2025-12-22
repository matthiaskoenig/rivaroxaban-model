from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import RivaroxabanSimulationExperiment
from pkdb_models.models.rivaroxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, \
    Fasting, RivaroxabanMappingMetaData, Coadministration
from pkdb_models.models.rivaroxaban.helpers import run_experiments

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim


class Moore2014(RivaroxabanSimulationExperiment):
    """Simulation experiment of Moore2014."""

    groups = ["NRF", "mildRI", "moderateRI"]

    group_labels = {
        "NRF": "Normal renal function",
        "mildRI": "Mild renal impairment",
        "moderateRI": "Moderate renal impairment",
    }

    colors = {
        "NRF": "black",
        "mildRI": "#66c2a4",
        "moderateRI": "#2ca25f",
    }

    bodyweights = {
        "NRF": 75.5,
        "mildRI": 76.7,
        "moderateRI": 84.8,
    }

    renal_functions = {
        "NRF": "Normal renal function",
        "mildRI": "Mild renal impairment",
        "moderateRI": "Moderate renal impairment",
    }

    infos = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Tab2A"]:
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
            rf_value = RivaroxabanSimulationExperiment.renal_map[self.renal_functions[group]]
            tcsims[group] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=50 * 60,
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweights[group], "kg"),
                        "PODOSE_riv": Q_(10, "mg"),
                        "KI__f_renal_function": Q_(rf_value, "dimensionless"),
                    },
                )
            ])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for kp, sid in enumerate(self.infos):
            name = self.infos[sid]
            for group in self.groups:
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
                        health=Health.HEALTHY if "NRF" in group else Health.RENAL_IMPAIRMENT,
                        fasting=Fasting.FED,
                        coadministration=Coadministration.NONE,
                    ),
                )
        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}
        fig = Figure(
            experiment=self,
            sid="Fig2_Tab2",
            name=self.__class__.__name__,
            num_rows=2,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)

        for kp, sid in enumerate(self.infos):
            name = self.infos[sid]
            for group in self.groups:
                plots[kp].add_data(
                    task=f"task_{group}",
                    xid="time",
                    yid=sid,
                    label = self.group_labels[group],
                    color=self.colors[group],
                )
                plots[kp].add_data(
                    dataset=f"{name}_{group}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label = self.group_labels[group],
                    linestyle="" if "urine" in name else "--",
                    color=self.colors[group],
                )

        figures[fig.sid] = fig

        return figures


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Moore2014.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Moore2014, output_dir=Moore2014.__name__)
