"""Simulation experiment of the study Kreutz2017."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Dosing,
    Fasting,
    Health,
    PKPDMappingMetaData,
    Route,
    Tissue,
)
from sbmlsim import Q
from sbmlsim.data import DataSet
from sbmlsim.fit import FitData, FitMapping
from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Simulation

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)


class Kreutz2017(RivaroxabanSimulationExperiment):
    """Simulation experiment of Kreutz2017.

    The rivaroxaban period of the crossover study: 20 mg once daily for 7 days
    with a standardized breakfast (fed) in healthy men. Plasma (Fig2, geometric
    mean and SD) and the prothrombin time and aPTT ratios to the baseline
    (Fig6) after the dose of day 7, the time is relative to that dose.
    """

    bodyweight = 81.1  # [kg]
    bodyheight = 181  # [cm]

    dose = 20  # [mg]
    route = "PO"

    info_figpk: ClassVar[dict[str, str]] = {
        "riv_rivaroxaban": "[Cve_riv]",
    }

    info_figpd: ClassVar[dict[str, str]] = {
        "riv_prothrombin time ratio": "PT_ratio",
        "riv_aPTT ratio": "aPTT_ratio",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig6"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if label.startswith("api_"):
                    continue
                dset = DataSet.from_df(df_label, self.ureg)
                if df_label["measurement_type"].iloc[0] == "concentration":
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                # the data are curated at the time after the first dose
                dset.time = dset.time - 6 * 24
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}

        # a dose every 24 hr for 7 days, the time is relative to the last dose;
        # a dose adds to what is left of the earlier ones
        tcsims["rivaroxaban"] = Simulation(
            time_unit="hr",
            end=6 * 24 + 73,
            preinit_changes={
                **self.default_changes(),
                "BW": Q(self.bodyweight, "kg"),
                f"{self.route}DOSE_riv": Q(self.dose, "mg"),
                "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
            },
            changes=[
                self.add_doses(
                    [24 * k for k in range(1, 7)],
                    {f"{self.route}DOSE_riv": Q(self.dose, "mg")},
                )
            ],
            time_shift=-6 * 24,
        )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}

        for variables in [self.info_figpk.items(), self.info_figpd.items()]:
            for label, sid in variables:
                mappings[self.to_sid(f"fm_{self.route}_{label}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(label),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task="task_rivaroxaban",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_pk(),
            **self._figure_pd(),
        }

    def _figure_pk(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="PK",
            name=self.__class__.__name__,
            num_rows=1,
            num_cols=1,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24, max=75),
            legend=True,
        )
        plots[0].set_yaxis(
            self.labels["[Cve_riv]"], unit=self.units["[Cve_riv]"], scale="linear"
        )

        for label, sid in self.info_figpk.items():
            # simulation
            plots[0].add_data(
                task="task_rivaroxaban",
                xid="time",
                yid=sid,
                label=f"{self.dose} mg QD, day 7 (fed)",
                color=self.fasting_colors["fed"],
            )
            # data
            plots[0].add_data(
                dataset=self.to_sid(label),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{self.dose} mg QD, day 7 (fed)",
                color=self.fasting_colors["fed"],
            )

        return self.figures_by_sid(fig)

    def _figure_pd(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="PD",
            name=self.__class__.__name__,
            num_rows=1,
            num_cols=2,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-24, max=75),
            legend=True,
        )

        for kp, (label, sid) in enumerate(self.info_figpd.items()):
            plots[kp].set_yaxis(self.labels[sid], unit=self.units[sid], scale="linear")
            # simulation
            plots[kp].add_data(
                task="task_rivaroxaban",
                xid="time",
                yid=sid,
                label=f"{self.dose} mg QD, day 7 (fed)",
                color=self.fasting_colors["fed"],
            )
            # data
            plots[kp].add_data(
                dataset=self.to_sid(label),
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=f"{self.dose} mg QD, day 7 (fed)",
                color=self.fasting_colors["fed"],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Kreutz2017, output_dir=Kreutz2017.__name__)
