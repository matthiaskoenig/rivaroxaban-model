from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import RivaroxabanSimulationExperiment
from pkdb_models.models.rivaroxaban.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, RivaroxabanMappingMetaData
)
from pkdb_models.models.rivaroxaban.helpers import run_experiments


class Kubitza2005a(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kubitza 2005a with mixed dosing (single + BID)."""

    formulations = ["T1.25", "T5", "T10", "T15", "T20", "T30", "T40", "T60", "T80", "S5", "S10"]
    pd_formulations = ["P", "T1.25", "T5", "T10", "T20", "T40", "T80"]
    bodyweight = 81.2

    colors = {
        "T1.25": "#FFFACD",  # Lemon Chiffon (very light yellow)
        "T5":     "#FFD580",  # Light orange (peachy)
        "T10":    "#FFB347",  # Soft orange
        "T15":    "#FF944D",  # Orange with a red tinge
        "T20":    "#FF7043",  # Reddish-orange
        "T30":    "#D25100",  # Burnt orange
        "T40":    "#A63C00",  # Darker brownish orange
        "T60":    "#7B2E00",  # Deep brown
        "T80":    "#4C1F00",  # Very dark brown

        "S5":     "#FFD580",
        "S10":    "#FFB347",

        "P":      "black",    # Placebo
    }

    formulation_labels = {
        "P": "Placebo",
        "T1.25": "Rivaroxaban 1.25 mg",
        "T5": "Rivaroxaban 5 mg",
        "T10": "Rivaroxaban 10 mg",
        "T15": "Rivaroxaban 15 mg",
        "T20": "Rivaroxaban 20 mg",
        "T30": "Rivaroxaban 30 mg",
        "T40": "Rivaroxaban 40 mg",
        "T60": "Rivaroxaban 60 mg",
        "T80": "Rivaroxaban 80 mg",
        "S5": "Rivaroxaban 5 mg Solution",
        "S10": "Rivaroxaban 10 mg Solution",
    }

    info_fig5 = {
        "[Cve_riv]": "rivaroxaban",
    }

    info_fig3 = {
        "Xa_inhibition": "factor Xa inhibition",
        "PT_ratio": "prothrombin time (change relative)",
        "aPTT_ratio": "aPTT (change relative)",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig5", "Fig3A", "Fig3B", "Fig3C"]

        for fig_id in fig_ids:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[f"{fig_id}_{label}"] = dset

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        all_formulations = list(set(self.formulations + self.pd_formulations))

        for f in all_formulations:
            dose = 0.0 if f == "P" else float(f[1:])
            tcsims[f] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=50 * 60,
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "PODOSE_riv": Q_(dose, "mg"),
                    }
                )
            ])

        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        for f in self.formulations:
            for sid, name in self.info_fig5.items():
                dataset_id = f"Fig5_{name}_{f}"
                mappings[f"fm_{dataset_id}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_id,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count"
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{f}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET if f.startswith("T") else ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED
                    )
                )

        for f in self.pd_formulations:
            for fig_id, (sid, name) in zip(["Fig3A", "Fig3B", "Fig3C"], self.info_fig3.items()):
                dataset_id = f"{fig_id}_{name}_{f}"
                mappings[f"fm_{dataset_id}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=dataset_id,
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count"
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{f}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=RivaroxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET if f.startswith("T") else ApplicationForm.SOLUTION,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED
                    )
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        figures = {}

        tablet_formulations = [f for f in self.formulations if f.startswith("T")]

        fig_pk1 = Figure(
            experiment=self, sid="Fig5",
            name=self.__class__.__name__
        )
        plots = fig_pk1.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)

        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        for f in tablet_formulations:
            plots[0].add_data(task=f"task_{f}", xid="time", yid="[Cve_riv]", label=self.formulation_labels[f], color=self.colors[f])
            plots[0].add_data(dataset=f"Fig5_rivaroxaban_{f}", xid="time", yid="mean", yid_sd="mean_sd", count="count", label=self.formulation_labels[f], color=self.colors[f])

        figures[fig_pk1.sid] = fig_pk1

        fig_pd = Figure(
            experiment=self, sid="Fig3", num_rows=3,
            name=self.__class__.__name__
        )
        plots = fig_pd.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)

        for i, (sid, label) in enumerate(self.info_fig3.items()):
            plots[i].set_yaxis(self.labels[sid], unit=self.units[sid])
            for f in self.pd_formulations:
                plots[i].add_data(task=f"task_{f}", xid="time", yid=sid, label=self.formulation_labels[f],
                                  color=self.colors[f])
                plots[i].add_data(dataset=f"Fig3{chr(65 + i)}_{label}_{f}", xid="time", yid="mean",
                                  count="count", label=self.formulation_labels[f], color=self.colors[f])

        figures[fig_pd.sid] = fig_pd


        # Tablet and Solution
        fig_pk2 = Figure(
            experiment=self, sid="FigX", num_rows=2,
            name=self.__class__.__name__
        )
        plots = fig_pk2.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)
        for kp in range(2):
            plots[kp].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[0].set_title("Solution (single dose)")
        for f in ["S5", "S10"]:
            plots[0].add_data(task=f"task_{f}", xid="time", yid="[Cve_riv]", label=self.formulation_labels[f], color=self.colors[f])
            plots[0].add_data(dataset=f"Fig5_rivaroxaban_{f}", xid="time", yid="mean", yid_sd="mean_sd", count="count", label=self.formulation_labels[f], color=self.colors[f])

        plots[1].set_title("Tablet (single dose)")
        for f in ["T5", "T10"]:
            plots[1].add_data(task=f"task_{f}", xid="time", yid="[Cve_riv]", label=self.formulation_labels[f] + " BID", color=self.colors[f])
            plots[1].add_data(dataset=f"Fig5_rivaroxaban_{f}", xid="time", yid="mean", yid_sd="mean_sd", count="count", label=self.formulation_labels[f] + " BID", color=self.colors[f])

        figures[fig_pk2.sid] = fig_pk2

        return figures


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Kubitza2005a.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Kubitza2005a, output_dir=Kubitza2005a.__name__)
