"""Simulation experiment of the study Frost2014."""

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


class Frost2014(RivaroxabanSimulationExperiment):
    """Simulation experiment of Frost2014 (the rivaroxaban arm).

    14 healthy subjects, rivaroxaban 10 mg once daily for 4 days in the clinic
    (and apixaban 2.5 mg twice daily in the other period of the crossover),
    food not reported. Fig2 is the arithmetic mean at steady state after the
    morning dose of day 4 without an error (n = 10 at 48 hr, n = 2 at 72 hr).
    The anti-FXa activity of Fig3 is in LMWH IU/ml, which the model does not
    describe.
    """

    bodyweight = 75.9  # [kg]
    #: the subjects of the study; the count of the data falls to 10 at 48 hr
    #: and 2 at 72 hr, the count of a mapping is the one of the study
    n_subjects = 14

    colors: ClassVar[dict[str, str]] = {"mRIV10": "black"}
    interventions: ClassVar[list[str]] = list(colors.keys())
    doses: ClassVar[dict[str, float]] = {"mRIV10": 10}
    infos_pk: ClassVar[dict[str, str]] = {"[Cve_riv]": "rivaroxaban"}

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for intervention in self.interventions:
            dose = self.doses[intervention]
            # a dose every 24 hr for 4 days, the time is relative to the last
            # dose; a dose adds to what is left of the earlier ones
            tcsims[intervention] = Simulation(
                time_unit="hr",
                end=3 * 24 + 75,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweight, "kg"),
                    "PODOSE_riv": Q(dose, "mg"),
                },
                changes=[self.add_doses([24, 48, 72], {"PODOSE_riv": Q(dose, "mg")})],
                time_shift=-3 * 24,
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            # pharmacokinetics
            for sid, name in self.infos_pk.items():
                mappings[f"fm_{name}_{intervention}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count=self.n_subjects,
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self._figure_pk(),
        }

    def _figure_pk(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time, min=-74, max=74),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=1e-3
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False

        for intervention in self.interventions:
            for sid, name in self.infos_pk.items():
                for plot in plots:
                    # simulation
                    plot.add_data(
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                        label="10 mg once daily",
                        color=self.colors[intervention],
                    )
                    # data
                    plot.add_data(
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count=self.n_subjects,
                        label="10 mg once daily",
                        color=self.colors[intervention],
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Frost2014, output_dir=Frost2014.__name__)
