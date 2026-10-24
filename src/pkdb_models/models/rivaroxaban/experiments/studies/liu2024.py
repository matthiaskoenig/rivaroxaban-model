"""Simulation experiment of the study Liu2024."""

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


class Liu2024(RivaroxabanSimulationExperiment):
    """Simulation experiment of Liu2024.

    Single doses of rivaroxaban 10, 15 and 20 mg (Xarelto tablets, the
    reference product of bioequivalence trials) in healthy Chinese subjects,
    fasted and fed; parallel groups, each subject once. The body weight is not
    reported (BMI 22-23 kg/m^2), the default of the model is used.
    """

    doses: ClassVar[list[int]] = [10, 15, 20]
    fasting_states: ClassVar[list[str]] = ["fasted", "fed"]

    fed: ClassVar[dict[str, float]] = {
        "fasted": RivaroxabanSimulationExperiment.fasting_map["fasted"],
        "fed": RivaroxabanSimulationExperiment.fasting_map["fed"],
    }

    colors: ClassVar[dict[str, str]] = {
        "fasted_10": "#DDB8EA",
        "fasted_15": "#DDB8EA",
        "fasted_20": "#DDB8EA",
        "fed_10": "#A35AC7",
        "fed_15": "#A35AC7",
        "fed_20": "#A35AC7",
    }

    info_figs: ClassVar[dict[str, dict[str, str]]] = {
        "Fig2": {"[Cve_riv]": "rivaroxaban"},
        "Fig3": {
            "PT_ratio": "prothrombin time (change relative)",
            "aPTT_ratio": "aPTT (change relative)",
        },
    }

    def cv_cmax(self) -> dict[str, float]:
        """Coefficient of variation of Cmax of every group (Tab 2)."""
        df = self.load_dataframe("Tab2")
        df = df[df.measurement_type == "cmax"]
        return {
            str(row.group): float(row.sd) / float(row.mean) for row in df.itertuples()
        }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        cv_cmax = self.cv_cmax()
        for fig_id in ["Fig2", "Fig3"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                if "Anti-factor Xa activity" in label:
                    # rivaroxaban by anti-Xa activity of the subjects of Fig2
                    continue
                if "digitised" in label:
                    continue  # take only table values
                if fig_id == "Fig2":
                    # the error bars of Fig2 (2.5 to 30 ng/ml) are not the SDs
                    # of Tab 2 (Cmax 44 to 76 ng/ml) and are smaller than them
                    # by 2 to 5 fold; the points get the coefficient of
                    # variation of Cmax of their group instead
                    group = label.removeprefix("rivaroxaban_")
                    df_label = df_label.copy()
                    df_label["sd"] = df_label["mean"] * cv_cmax[group]
                    df_label["se"] = float("nan")
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(f"{fig_id}_{label}")] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for fasting in self.fasting_states:
            for dose in self.doses:
                sim_id = f"{fasting}_{dose}"
                simulations[self.to_sid(sim_id)] = Simulation(
                    time_unit="min",
                    end=50 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "PODOSE_riv": Q(dose, "mg"),
                        "GU__fed": Q(self.fed[fasting], "dimensionless"),
                    },
                )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for fig_id, variables in self.info_figs.items():
            for yid, label in variables.items():
                for fasting in self.fasting_states:
                    for dose in self.doses:
                        sim_id = f"{fasting}_{dose}"
                        dataset_id = f"{fig_id}_{label}_{sim_id}"

                        mappings[self.to_sid(f"fm_{dataset_id}")] = FitMapping(
                            self,
                            reference=FitData(
                                self,
                                dataset=self.to_sid(dataset_id),
                                xid="time",
                                yid="mean",
                                yid_sd="mean_sd",
                                count="count",
                            ),
                            observable=FitData(
                                self,
                                task=self.to_sid(f"task_{sim_id}"),
                                xid="time",
                                yid=yid,
                            ),
                            metadata=PKPDMappingMetaData(
                                tissue=Tissue.PLASMA,
                                route=Route.PO,
                                application_form=ApplicationForm.TABLET,
                                dosing=Dosing.SINGLE,
                                health=Health.HEALTHY,
                                fasting=Fasting.FASTED
                                if fasting == "fasted"
                                else Fasting.FED,
                                coadministration=Coadministration.NONE,
                            ),
                        )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        return {
            **self.figures_pk_pd(),
        }

    def figures_pk_pd(self) -> dict[str, Figure]:
        """Figures of the plasma concentrations and the coagulation tests."""
        # PK Figure: Fig2
        fig_pk = Figure(
            experiment=self,
            sid="Fig2",
            name=self.__class__.__name__,
            num_cols=3,
        )
        plots_pk = fig_pk.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        for kp in range(3):
            plots_pk[kp].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)

        for fasting in self.fasting_states:
            for kp, dose in enumerate(self.doses):
                sim_id = f"{fasting}_{dose}"
                plots_pk[kp].add_data(
                    task=self.to_sid(f"task_{sim_id}"),
                    xid="time",
                    yid="[Cve_riv]",
                    label=f"RIV {dose} mg {fasting}",
                    color=self.colors[sim_id],
                )
                plots_pk[kp].add_data(
                    dataset=self.to_sid(f"Fig2_rivaroxaban_{sim_id}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"RIV {dose} mg {fasting}",
                    color=self.colors[sim_id],
                )

        # PD Figure: Fig3 (combined PT_ratio and aPTT_ratio)
        fig_pd = Figure(
            experiment=self,
            sid="Fig3",
            num_rows=2,
            num_cols=3,
            name=self.__class__.__name__,
        )
        plots_pd = fig_pd.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )

        for kp in range(3):
            plots_pd[kp].set_yaxis(self.labels["PT_ratio"], unit=self.units["PT_ratio"])
            plots_pd[kp + 3].set_yaxis(
                self.labels["aPTT_ratio"], unit=self.units["aPTT_ratio"]
            )

        for fasting in self.fasting_states:
            for kp, dose in enumerate(self.doses):
                sim_id = f"{fasting}_{dose}"
                color = self.colors[sim_id]

                # Subplot 1: PT_ratio
                plots_pd[kp].add_data(
                    task=self.to_sid(f"task_{sim_id}"),
                    xid="time",
                    yid="PT_ratio",
                    label=f"RIV {dose} mg {fasting}",
                    color=color,
                )
                plots_pd[kp].add_data(
                    dataset=self.to_sid(
                        f"Fig3_prothrombin time (change relative)_{sim_id}"
                    ),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"RIV {dose} mg {fasting}",
                    color=color,
                )

                # Subplot 2: aPTT_ratio
                plots_pd[kp + 3].add_data(
                    task=self.to_sid(f"task_{sim_id}"),
                    xid="time",
                    yid="aPTT_ratio",
                    label=f"RIV {dose} mg {fasting}",
                    color=color,
                )
                plots_pd[kp + 3].add_data(
                    dataset=self.to_sid(f"Fig3_aPTT (change relative)_{sim_id}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"RIV {dose} mg {fasting}",
                    color=color,
                )

        return self.figures_by_sid(fig_pk, fig_pd)


if __name__ == "__main__":
    MODEL.run_experiments(Liu2024, output_dir=Liu2024.__name__)
