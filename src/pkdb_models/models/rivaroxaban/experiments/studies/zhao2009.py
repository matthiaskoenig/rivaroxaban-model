from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim
from sbmlutils.console import console

from pkdb_models.models import rivaroxaban
from pkdb_models.models.rivaroxaban.experiments.base_experiment import RivaroxabanSimulationExperiment
from pkdb_models.models.rivaroxaban.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health, Fasting, RivaroxabanMappingMetaData, Coadministration
)
from pkdb_models.models.rivaroxaban.helpers import run_experiments


class Zhao2009(RivaroxabanSimulationExperiment):
    """Simulation experiment of Zhao 2009."""

    formulations = [
        "SDP", "SD2.5", "SD5", "SD10", "SD20", "SD40",
        "MDP", "MD5", "MD10", "MD20", "MD30"
    ]

    bodyweights = {
        "SDP": 63.8, "SD2.5": 62.4, "SD5": 66, "SD10": 58.4, "SD20": 59.1, "SD30": 62.1, "SD40": 62.8,
        "MDP": 66.3, "MD5": 66, "MD10": 64, "MD20": 68.3, "MD30": 65.6,
    }

    colors = {
        # All doses use a unified gradient regardless of SD/MD
        "SDP": "black",
        "SD2.5": "#FFFACD",  # Very light yellow (like T1.25)
        "SD5": "#FFE299",  # Pale warm yellow
        "SD10": "#FFB347",  # Soft orange (like T10)
        "SD20": "#FF944D",  # Reddish orange (like T15–T20)
        "SD30": "#D2691E",  # Chocolate/burnt orange (like T30)
        "SD40": "#8B4513",  # Saddle brown (like T40–T60)

        "MDP": "black",
        "MD5": "#FFE299",  # Same as SD5
        "MD10": "#FFB347",  # Same as SD10
        "MD20": "#FF944D",  # Same as SD20
        "MD30": "#D2691E",  # Same as SD30
    }

    info_fig_pk = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine"
    }

    info_fig_pd = {
        "Xa_inhibition": "factor Xa inhibition",
        "PT_ratio": "prothrombin time (change relative)",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        fig_ids = ["Fig2C", "Fig1A_2A", "Fig1B_2B", "Tab2_4U"]

        for fig_id in fig_ids:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[f"{label}"] = dset

        # console.print(dsets.keys())
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        for formulation in self.formulations:
            dose = 0 if formulation.endswith("P") else int(float(formulation[2:]))
            bw = self.bodyweights.get(formulation)

            if formulation.startswith("SD"):
                # Single dose
                tcsims[formulation] = TimecourseSim([
                    Timecourse(
                        start=0,
                        end=52 * 60,
                        steps=1000,
                        changes={
                            **self.default_changes(),
                            "BW": Q_(bw, "kg"),
                            "PODOSE_riv": Q_(dose, "mg"),
                        }
                    )
                ])

            elif formulation.startswith("MD"):
                # Multiple dose BID for 6 days

                tc0 = Timecourse(
                    start=0,
                    end=48 * 60,
                    steps=300,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(bw, "kg"),
                        "PODOSE_riv": Q_(dose, "mg")
                    },
                )
                tc1 = Timecourse(
                    start=0,
                    end=12 * 60,
                    steps=300,
                    changes={
                        "PODOSE_riv": Q_(dose, "mg"),
                        "Aurine_riv": Q_(0, "mmole"),
                    },
                )
                tc2 = Timecourse(
                    start=0,
                    end=52 * 60,
                    steps=300,
                    changes={
                        "PODOSE_riv": Q_(dose, "mg"),
                        "Aurine_riv": Q_(0, "mmole"),
                    },
                )

                tcsims[formulation] = TimecourseSim(
                    [tc0] + [tc1 for _ in range(6)] + [tc2],
                    # time_offset=-5*24*60,
                )


        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}

        # pharmacokinetics
        for krow, sid in enumerate(self.info_fig_pk):
            name = self.info_fig_pk[sid]
            for formulation in self.formulations:
                if formulation.endswith("P"):
                    # no placebo PK
                    continue

                if formulation.startswith("SD") and not "urine" in name:
                    # SD has only urine PK
                    continue

                mappings[f"fm_{name}_{formulation}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{formulation}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if "urine" in name else None,
                        count="count"
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
                        dosing=Dosing.SINGLE if formulation.startswith("SD") else Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    )
                )

        # pharmacodynamics
        for krow, sid in enumerate(self.info_fig_pd):
            name = self.info_fig_pd[sid]
            for formulation in self.formulations:
                mappings[f"fm_{name}_{formulation}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{formulation}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count"
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
                        dosing=Dosing.SINGLE if formulation.startswith("SD") else Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE,
                    )
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        fig = Figure(experiment=self, sid="Fig2C", num_rows=1, num_cols=3, name=self.__class__.__name__)
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)

        plots[0].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
        plots[1].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[2].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)

        for krow, sid in enumerate(self.info_fig_pk):
            name = self.info_fig_pk[sid]
            for formulation in self.formulations:

                if formulation.startswith("SD"):
                    if "urine" in name:
                        kp = 0
                    else:
                        continue
                elif formulation.startswith("MD"):
                    kp = 2 if "urine" in name else 1

                # simulation
                plots[kp].add_data(
                    task=f"task_{formulation}",
                    xid="time",
                    yid=sid,
                    label=formulation,
                    color=self.colors[formulation]
                )

                # data
                if formulation.endswith("P"):
                    continue
                plots[kp].add_data(
                    dataset=f"{name}_{formulation}",
                    xid="time", yid="mean",
                    yid_sd="mean_sd" if "urine" in name else None,
                    count="count",
                    label=formulation,
                    linestyle="" if "urine" in name else "--",
                    color=self.colors[formulation]
                )

        return {fig.sid: fig}

    def figure_pd(self) -> Dict[str, Figure]:
        fig = Figure(experiment=self, sid="Fig1", num_rows=2, num_cols=2, name=self.__class__.__name__)
        plots = fig.create_plots(xaxis=Axis(self.label_time, unit=self.unit_time), legend=True)

        # Xa inhibition
        plots[0].set_yaxis(self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"])
        # Xa inhibition - multiple dose
        plots[1].set_yaxis(self.labels["Xa_inhibition"], unit=self.units["Xa_inhibition"])

        # PT ratio
        plots[2].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
        # PT ratio - multiple dose
        plots[3].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])

        for formulation in self.formulations:
            is_sd = formulation.startswith("SD") or formulation == "SDP"
            col_idx = 0 if is_sd else 1  # Xa_inhibition plots
            row_idx = 0  # First row

            plots[row_idx * 2 + col_idx].add_data(
                task=f"task_{formulation}",
                xid="time",
                yid="Xa_inhibition",
                label=formulation,
                color=self.colors.get(formulation, "black")
            )
            plots[row_idx * 2 + col_idx].add_data(
                dataset=f"factor Xa inhibition_{formulation}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=formulation,
                color=self.colors.get(formulation, "black")
            )

            col_idx = 0 if is_sd else 1  # PT_ratio plots
            row_idx = 1  # Second row

            plots[row_idx * 2 + col_idx].add_data(
                task=f"task_{formulation}",
                xid="time",
                yid="PT_ratio",
                label=formulation,
                color=self.colors.get(formulation, "black")
            )
            plots[row_idx * 2 + col_idx].add_data(
                dataset=f"prothrombin time (change relative)_{formulation}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=formulation,
                color=self.colors.get(formulation, "black")
            )

        return {fig.sid: fig}


if __name__ == "__main__":
    out = rivaroxaban.RESULTS_PATH_SIMULATION / Zhao2009.__name__
    out.mkdir(parents=True, exist_ok=True)
    run_experiments(Zhao2009, output_dir=Zhao2009.__name__)