"""Simulation experiment of the study Metze2021."""

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


class Metze2021(RivaroxabanSimulationExperiment):
    """Simulation experiment of Metze2021.

    20 patients (age 72 +- 13 years, atrial fibrillation or venous
    thromboembolism) at the start of rivaroxaban 20 mg once daily: the plasma
    concentrations (anti-Xa activity with rivaroxaban calibrators) before and
    3, 6 and 12 hr after the first dose of a course. The mean creatinine
    clearance of 67 +- 18 ml/min is simulated as a mild renal impairment
    (67/101); the body weight is not reported.
    """

    color = "#66c2a4"  # Mild renal impairment
    #: the 20 individuals IS041 to IS060 (Fig1, vector paths of the PDF)
    individuals: ClassVar[list[str]] = [f"IS{k:03d}" for k in range(41, 61)]
    #: creatinine clearance of the rivaroxaban group [ml/min] (Tab 1)
    crcl = 67.0
    intervention = "RIV20"
    dose = 20  # [mg]

    infos_pk: ClassVar[dict[str, str]] = {"[Cve_riv]": "rivaroxaban"}

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig1A"]:
            df = self.load_dataframe(fig_id)
            df = df[df.substance == "rivaroxaban"]
            if fig_id == "Fig1A":
                # the levels before the dose (up to 55 ng/ml) are left from
                # earlier intakes, which the simulation of a first dose has not
                df = df[df.time > 0]
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion(
                    "value" if fig_id == "Fig1" else "mean", 1 / self.Mr.riv
                )
                dsets[label] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        simulations[self.intervention] = Simulation(
            time_unit="min",
            end=16 * 60,
            preinit_changes={
                **self.default_changes(),
                "PODOSE_riv": Q(self.dose, "mg"),
                # renal function: the creatinine clearance relative to the
                # GFR of 101 ml/min of the model
                "KI__f_renal_function": Q(self.crcl / 101, "dimensionless"),
            },
        )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.infos_pk.items():
            # PK (individual)
            for individual in self.individuals:
                mappings[f"fm_{name}_{self.intervention}_{individual}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{self.intervention}_{individual}",
                        xid="time",
                        yid="value",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{self.intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.RENAL_IMPAIRMENT,
                        fasting=Fasting.NR,
                        coadministration=Coadministration.NONE,
                    ),
                )

            # PK (mean)
            mappings[f"fm_{name}_{self.intervention}"] = FitMapping(
                self,
                reference=FitData(
                    self,
                    dataset=f"{name}_{self.intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                ),
                observable=FitData(
                    self,
                    task=f"task_{self.intervention}",
                    xid="time",
                    yid=sid,
                ),
                metadata=PKPDMappingMetaData(
                    tissue=Tissue.PLASMA,
                    route=Route.PO,
                    application_form=ApplicationForm.TABLET,
                    dosing=Dosing.SINGLE,
                    health=Health.RENAL_IMPAIRMENT,
                    fasting=Fasting.NR,
                    coadministration=Coadministration.NONE,
                ),
            )

        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> dict[str, Figure]:
        """Figure of the individual and mean plasma concentrations."""
        fig = Figure(
            experiment=self, sid="Fig1", name=f"{self.__class__.__name__} (patients)"
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for ks, (sid, name) in enumerate(self.infos_pk.items()):
            # individual data
            for ki, individual in enumerate(self.individuals):
                plots[ks].add_data(
                    dataset=f"{name}_{self.intervention}_{individual}",
                    xid="time",
                    yid="value",
                    yid_sd=None,
                    count="count",
                    label="RIV 20 mg, individuals" if ki == 0 else None,
                    color=self.color,
                )

            # simulation
            plots[ks].add_data(
                task=f"task_{self.intervention}",
                xid="time",
                yid=sid,
                label="RIV 20 mg",
                color="black",
            )

            # mean data
            plots[ks].add_data(
                dataset=f"{name}_{self.intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label="RIV 20 mg",
                color="black",
                linestyle="-",
            )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Metze2021, output_dir=Metze2021.__name__)
