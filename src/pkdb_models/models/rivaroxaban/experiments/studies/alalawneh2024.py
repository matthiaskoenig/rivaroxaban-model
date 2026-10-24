"""Simulation experiment of the study Alalawneh2024."""

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


class Alalawneh2024(RivaroxabanSimulationExperiment):
    """Simulation experiment of Alalawneh2024 (RIVOBESE-PK).

    36 healthy men, 18 obese (BMI >= 35) and 18 non-obese, a single 20 mg
    tablet of Xarelto 30 min after the start of a high-fat breakfast (fed).
    The morbidly obese (BMI > 40 or BW > 120 kg, n = 6, Fig3) are a subset of
    the obese group. PT and aPTT (Fig4) are absolute times of the laboratory
    without an error, the baselines of the model are set to the predose values
    of each group.
    """

    groups: ClassVar[list[str]] = ["non-obese", "obese", "morbidly-obese"]
    bodyweights: ClassVar[dict[str, float]] = {
        "non-obese": 67.4,  # kg
        "obese": 111.5,  # kg
        # the heaviest 6 of the 18 obese (111.5 +- 11.7 kg): the expected mean
        # of the top third of a normal distribution, 111.5 + 1.09 * 11.7 kg
        "morbidly-obese": 124,  # kg
    }
    health: ClassVar[dict[str, Health]] = {
        "non-obese": Health.HEALTHY,
        "obese": Health.OBESE,
        "morbidly-obese": Health.OBESE,
    }
    #: predose PT and aPTT of the groups (Fig4a/b), the morbidly obese take
    #: the values of the obese group
    pt_ref: ClassVar[dict[str, float]] = {
        "non-obese": 12.93,  # s
        "obese": 12.77,  # s
        "morbidly-obese": 12.77,  # s
    }
    aptt_ref: ClassVar[dict[str, float]] = {
        "non-obese": 40.32,  # s
        "obese": 35.13,  # s
        "morbidly-obese": 35.13,  # s
    }
    colors: ClassVar[dict[str, str]] = {
        "non-obese": "black",
        "obese": "tab:blue",
        "morbidly-obese": "tab:orange",
    }
    info_fig3: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
        "KI__RIVEX": "rivaroxaban_excretion_urine",
    }
    info_fig4: ClassVar[dict[str, str]] = {
        "PT": "prothrombin time",
        "aPTT": "aPTT",
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3", "Fig4a", "Fig4b", "Tab4"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                if "rivaroxaban" in label:
                    dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[self.to_sid(label)] = dset

        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        tcsims = {}
        for group in self.groups:
            tcsims[self.to_sid(f"{group}")] = Simulation(
                time_unit="min",
                end=50 * 60,
                preinit_changes={
                    **self.default_changes(),
                    "BW": Q(self.bodyweights[group], "kg"),
                    "PODOSE_riv": Q(20, "mg"),
                    "GU__fed": Q(self.fasting_map["fed"], "dimensionless"),
                    "PT_ref": Q(self.pt_ref[group], "s"),
                    "aPTT_ref": Q(self.aptt_ref[group], "s"),
                },
            )

        return tcsims

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            # pharmacokinetics
            for k, (sid, name) in enumerate(self.info_fig3.items()):
                if k != 0 and "morbid" in group:
                    continue
                mappings[self.to_sid(f"fm_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{group}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=self.health[group],
                        fasting=Fasting.FED,
                        coadministration=Coadministration.NONE,
                    ),
                )

            # pharmacodynamics
            if "morbid" in group:
                continue
            for sid, name in self.info_fig4.items():
                mappings[self.to_sid(f"fm_{name}_{group}")] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=self.to_sid(f"task_{group}"),
                        xid="time",
                        yid=sid,
                    ),
                    metadata=PKPDMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=self.health[group],
                        fasting=Fasting.FED,
                        coadministration=Coadministration.NONE,
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
            sid="Fig2_3_Tab4",
            num_rows=2,
            num_cols=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
        plots[1].set_yaxis(
            self.label_riv_plasma, unit=self.unit_riv, scale="log", min=1e-3
        )
        # the log panel shares the legend of the linear one
        plots[1].legend = False
        plots[2].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
        plots[3].set_yaxis(self.labels["KI__RIVEX"], unit=self.units["KI__RIVEX"])
        # the panels of the observables, the plasma is shown linear and log
        panels: dict[str, list[int]] = {
            "[Cve_riv]": [0, 1],
            "Aurine_riv": [2],
            "KI__RIVEX": [3],
        }

        for group in self.groups:
            for k, (sid, name) in enumerate(self.info_fig3.items()):
                if k != 0 and "morbid" in group:
                    continue
                for ip in panels[sid]:
                    plots[ip].add_data(
                        task=self.to_sid(f"task_{group}"),
                        xid="time",
                        yid=sid,
                        label=group,
                        color=self.colors[group],
                    )
                    plots[ip].add_data(
                        dataset=self.to_sid(f"{name}_{group}"),
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=group,
                        color=self.colors[group],
                    )

        return self.figures_by_sid(fig)

    def _figure_pd(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig4",
            num_rows=2,
            name=self.__class__.__name__,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.labels["PT"], unit=self.units["PT"])
        plots[1].set_yaxis(self.labels["aPTT"], unit=self.units["aPTT"])
        for group in self.groups:
            if "morbid" in group:
                continue
            for k, (sid, name) in enumerate(self.info_fig4.items()):
                plots[k].add_data(
                    task=self.to_sid(f"task_{group}"),
                    xid="time",
                    yid=sid,
                    label=group,
                    color=self.colors[group],
                )
                plots[k].add_data(
                    dataset=self.to_sid(f"{name}_{group}"),
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=group,
                    color=self.colors[group],
                )
        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Alalawneh2024, output_dir=Alalawneh2024.__name__)
