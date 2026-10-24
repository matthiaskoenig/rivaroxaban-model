"""Simulation experiment of the study GenisNajera2022."""

from typing import ClassVar, override

from pkpdsim.metadata import (
    ApplicationForm,
    Coadministration,
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


class GenisNajera2022(RivaroxabanSimulationExperiment):
    """Simulation experiment of GenisNajera2022.

    Bioequivalence of a generic (T) and Xarelto (R) 20 mg tablet in 31 healthy
    Mexican volunteers (29 female), a crossover with a light standardized
    breakfast (about 380 kcal, low fat) 0.5 hr before the dose, simulated as
    the fed state. T and R are the same subjects and are simulated alike.
    Fig1 is the mean +- SE.
    """

    formulations: ClassVar[list[str]] = ["T20", "R20"]
    colors: ClassVar[dict[str, str]] = {
        "R20": "black",
        "T20": "tab:blue",
    }
    bodyweight = 61.7  # [kg]

    formulation_labels: ClassVar[dict[str, str]] = {
        "T20": "test 20 mg",
        "R20": "Xarelto 20 mg",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig1")
        for label, df_label in df.groupby("label"):
            label = str(label)
            dset = DataSet.from_df(df_label, self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.riv)
            dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for formulation in self.formulations:
            tcsims[f"po_{formulation}"] = Simulation(
                time_unit="min",
                end=50 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "PODOSE_riv": Q(20, "mg"),
                    "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
                },
            )
        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
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
                    self,
                    task=f"task_po_{formulation}",
                    xid="time",
                    yid="[Cve_riv]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.LIGHT_BREAKFAST,
                    coadministration=Coadministration.NONE,
                ),
            )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=1e-3
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False

        for plot in plots:
            for formulation in self.formulations:
                # simulation
                plot.add_data(
                    task=f"task_po_{formulation}",
                    xid="time",
                    yid="[Cve_riv]",
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )

                # data
                plot.add_data(
                    dataset=f"rivaroxaban_{formulation}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=self.formulation_labels[formulation],
                    color=self.colors[formulation],
                )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(GenisNajera2022, output_dir=GenisNajera2022.__name__)
