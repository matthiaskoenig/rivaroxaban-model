"""Simulation experiment of the study Moore2014."""

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


class Moore2014(RivaroxabanSimulationExperiment):
    """Simulation experiment of Moore2014.

    Single dose of rivaroxaban 10 mg in subjects with normal renal function
    (n=8) and with mild (n=12) or moderate (n=9) renal impairment, alone and
    on day 5 of erythromycin 500 mg three times daily (n=8 per group), plasma
    and urine. Food is not reported.
    """

    groups: ClassVar[list[str]] = ["NRF", "mildRI", "moderateRI"]

    group_labels: ClassVar[dict[str, str]] = {
        "NRF": "normal renal function",
        "mildRI": "mild renal impairment",
        "moderateRI": "moderate renal impairment",
    }

    colors: ClassVar[dict[str, str]] = {
        "NRF": "black",
        "mildRI": "#66c2a4",
        "moderateRI": "#2ca25f",
    }
    color_ery = "tab:red"

    bodyweights: ClassVar[dict[str, float]] = {
        "NRF": 75.5,
        "mildRI": 76.7,
        "moderateRI": 84.8,
    }

    #: creatinine clearance (Cockcroft-Gault, mean of two) [ml/min] (Tab 1);
    #: the renal function is the clearance relative to the GFR of 101 ml/min
    crcl: ClassVar[dict[str, float]] = {
        "NRF": 99.5,
        "mildRI": 67.1,
        "moderateRI": 42.5,
    }

    #: rivaroxaban alone and with erythromycin
    arms: ClassVar[list[str]] = ["", "_ERY"]

    #: inhibition by erythromycin 500 mg three times daily (not fitted): a
    #: moderate CYP3A4 inhibitor (midazolam AUC about 3 to 4 fold, about 70 %
    #: of the CYP3A4 activity inhibited) and an inhibitor of the renal P-gp as
    #: clarithromycin (Lenard2024). Observed AUC ratios of rivaroxaban 1.38
    #: (normal), 1.60 (mild) and 1.53 (moderate renal impairment) in this
    #: study, 1.34 in Mueck2013.
    inhibition: ClassVar[dict[str, dict[str, float]]] = {
        "": {"LI__f_cyp3a4": 1.0, "KI__f_pgp": 1.0},
        "_ERY": {"LI__f_cyp3a4": 0.3, "KI__f_pgp": 0.7},
    }

    infos: ClassVar[dict[str, str]] = {
        "[Cve_riv]": "rivaroxaban",
        "Aurine_riv": "rivaroxaban_urine",
    }

    tissues: ClassVar[dict[str, Tissue]] = {
        "[Cve_riv]": Tissue.PLASMA,
        "Aurine_riv": Tissue.URINE,
    }

    @override
    def datasets(self) -> dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Tab2A"]:
            df = self.load_dataframe(fig_id)
            for label, df_label in df.groupby("label"):
                label = str(label)
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("mean", 1 / self.Mr.riv)
                dsets[label] = dset
        return dsets

    @override
    def simulations(self) -> dict[str, Simulation]:
        simulations = {}
        for group in self.groups:
            for arm in self.arms:
                simulations[f"{group}{arm}"] = Simulation(
                    time_unit="min",
                    end=50 * 60,
                    preinit_changes={
                        **self.default_changes(),
                        "BW": Q(self.bodyweights[group], "kg"),
                        "PODOSE_riv": Q(10, "mg"),
                        "KI__f_renal_function": Q(
                            self.crcl[group] / 101, "dimensionless"
                        ),
                        **{
                            pid: Q(value, "dimensionless")
                            for pid, value in self.inhibition[arm].items()
                        },
                    },
                )
        return simulations

    @override
    def fit_mappings(self) -> dict[str, FitMapping]:
        mappings = {}
        for sid, name in self.infos.items():
            for group in self.groups:
                for arm in self.arms:
                    mappings[f"fm_{name}_{group}{arm}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{group}{arm}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self,
                            task=f"task_{group}{arm}",
                            xid="time",
                            yid=sid,
                        ),
                        metadata=PKPDMappingMetaData(
                            tissue=self.tissues[sid],
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY
                            if group == "NRF"
                            else Health.RENAL_IMPAIRMENT,
                            fasting=Fasting.NR,
                            coadministration=Coadministration.ERYTHROMYCIN
                            if arm
                            else Coadministration.NONE,
                        ),
                    )
        return mappings

    @override
    def figures(self) -> dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2_Tab2",
            name=self.__class__.__name__,
            num_rows=2,
            num_cols=3,
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        for kg, group in enumerate(self.groups):
            plots[kg].set_title(self.group_labels[group])
            plots[kg].set_yaxis(self.label_riv_plasma, unit=self.unit_riv)
            plots[kg + 3].set_yaxis(self.label_riv_urine, unit=self.unit_riv_urine)
            for kp, (sid, name) in enumerate(self.infos.items()):
                plot = plots[kg + 3 * kp]
                for arm in self.arms:
                    label = "RIV 10 mg + erythromycin" if arm else "RIV 10 mg"
                    color = self.color_ery if arm else self.colors[group]
                    plot.add_data(
                        task=f"task_{group}{arm}",
                        xid="time",
                        yid=sid,
                        label=label,
                        color=color,
                    )
                    plot.add_data(
                        dataset=f"{name}_{group}{arm}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=label,
                        linestyle="" if "urine" in name else "--",
                        color=color,
                    )

        return self.figures_by_sid(fig)


if __name__ == "__main__":
    MODEL.run_experiments(Moore2014, output_dir=Moore2014.__name__)
