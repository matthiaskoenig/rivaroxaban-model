"""Simulation experiment of the study Lenard2024."""

from typing import ClassVar, override

import pandas as pd
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


def impute_sd(df: pd.DataFrame) -> pd.DataFrame:
    """Give the points of a mean curve without a standard deviation one.

    The figures of the microdose studies have no error bars after 8 hr (and
    hidden ones near the peak); sbmlsim gives a point without a standard
    deviation the largest one of its curve, which removes the low
    concentrations of the terminal phase from a fit. Such a point (after the
    dose) gets the median coefficient of variation of the points of its curve
    with a standard deviation.

    Args:
        df: the rows of one mean curve with the columns `time`, `mean`, `sd`.

    Returns:
        A copy of the rows with the standard deviations.
    """
    df = df.copy()
    has_sd = df["sd"].notna() & (df["mean"] > 0)
    cv = (df.loc[has_sd, "sd"] / df.loc[has_sd, "mean"]).median()
    missing = df["sd"].isna() & (df["time"] > 0)
    df.loc[missing, "sd"] = df.loc[missing, "mean"] * cv
    return df


class Lenard2024(RivaroxabanSimulationExperiment):
    """Simulation experiment of Lenard2024.

    Microdose cocktail of apixaban 25 µg, rivaroxaban 25 µg and edoxaban 50 µg
    as an oral solution in 12 healthy subjects, given alone and on day 8 of
    clarithromycin 500 mg twice daily (inhibition of CYP3A4 and P-gp at steady
    state). Only the rivaroxaban data is used.
    """

    interventions: ClassVar[list[str]] = [
        "API25, RIV25, EDO50",
        "API25, RIV25, EDO50, CLAR",
    ]
    #: keys of the datasets and mappings by intervention, the simulation is
    #: the key in upper case
    keys: ClassVar[dict[str, str]] = {
        "API25, RIV25, EDO50": "riv25",
        "API25, RIV25, EDO50, CLAR": "riv25_clar",
    }
    colors: ClassVar[dict[str, str]] = {
        "API25, RIV25, EDO50": "black",
        "API25, RIV25, EDO50, CLAR": "tab:red",
    }
    intervention_labels: ClassVar[dict[str, str]] = {
        "API25, RIV25, EDO50": "RIV 25 µg",
        "API25, RIV25, EDO50, CLAR": "RIV 25 µg + clarithromycin",
    }
    coadministrations: ClassVar[dict[str, Coadministration]] = {
        "API25, RIV25, EDO50": Coadministration.NONE,
        "API25, RIV25, EDO50, CLAR": Coadministration.CLARITHROMYCIN,
    }
    #: inhibition by clarithromycin 500 mg twice daily (not fitted): a strong
    #: CYP3A4 inhibitor (midazolam AUC about 8 fold, i.e. about 90 % of the
    #: CYP3A4 activity inhibited) and an inhibitor of the renal P-gp (renal
    #: clearance of digoxin reduced by about a third, Rengelshausen2003).
    #: Observed AUC ratio of rivaroxaban 1.43 (Cmax 1.40) in this study, 1.54
    #: with 20 mg (Mueck2013); the model has no intestinal CYP3A4 or P-gp, so
    #: the increase of Cmax is not described.
    inhibition: ClassVar[dict[str, dict[str, float]]] = {
        "API25, RIV25, EDO50": {"LI__f_cyp3a4": 1.0, "KI__f_pgp": 1.0},
        "API25, RIV25, EDO50, CLAR": {"LI__f_cyp3a4": 0.1, "KI__f_pgp": 0.7},
    }
    #: rivaroxaban 25 µg, a single dose in each period (the second dose on
    #: day 8 after a washout of 7 days is the clarithromycin arm)
    dose = 0.025  # [mg]

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        df = self.load_dataframe("Fig2")
        df = df[df.substance == "rivaroxaban"]
        for intervention, df_label in df.groupby("intervention"):
            dset = DataSet.from_df(impute_sd(df_label), self.ureg)
            dset.unit_conversion("mean", 1 / self.Mr.riv)
            dsets[self.keys[str(intervention)]] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        return {
            self.keys[intervention].upper(): Simulation(
                time_unit="min",
                end=25 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "POSOLUTION_riv": Q(self.dose, "mg"),
                    **{
                        pid: Q(value, "dimensionless")
                        for pid, value in self.inhibition[intervention].items()
                    },
                },
            )
            for intervention in self.interventions
        }

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            key = self.keys[intervention]
            mappings[f"fm_{key}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=key,
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_{key.upper()}",
                    xid="time",
                    yid="[Cve_riv]",
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.SOLUTION,
                    dosing=Dosing.SINGLE,
                    health=Health.HEALTHY,
                    fasting=Fasting.NR,
                    coadministration=self.coadministrations[intervention],
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=1,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for intervention in self.interventions:
            key = self.keys[intervention]
            plots[0].add_data(
                task=f"task_{key.upper()}",
                xid="time",
                yid="[Cve_riv]",
                label=self.intervention_labels[intervention],
                color=self.colors[intervention],
            )
            plots[0].add_data(
                dataset=key,
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.intervention_labels[intervention],
                color=self.colors[intervention],
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Lenard2024, output_dir=Lenard2024.__name__)
