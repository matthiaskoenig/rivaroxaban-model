"""Base class of the simulation experiments of the coagulation model alone."""

from collections.abc import Mapping
from pathlib import Path
from typing import ClassVar, override

from pkpdsim.experiment import PKPDSimulationExperiment
from sbmlsim.data import Data
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task
from sbmlsim.units import Quantity

from pkdb_models.models.rivaroxaban import MODEL, MODEL_BASE_PATH


class CoagulationSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the simulation experiments of the coagulation model alone.

    The plasma concentration of rivaroxaban is set by the simulations.
    """

    pkpd_model = MODEL

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 40

    units: ClassVar[Mapping[str, str]] = {
        "time": "hr",
        "[Cve_riv]": "µM",
        "PT": "s",
        "PT_change": "s",
        "PT_ratio": "dimensionless",
        "aPTT": "s",
        "aPTT_change": "s",
        "aPTT_ratio": "dimensionless",
        "Xa_inhibition": "dimensionless",
    }
    labels: ClassVar[Mapping[str, str]] = {
        "time": "time",
        "[Cve_riv]": "rivaroxaban plasma",
        "PT": "prothrombin time",
        "PT_change": "prothrombin time change",
        "PT_ratio": "prothrombin time ratio",
        "aPTT": "aPTT",
        "aPTT_change": "aPTT change",
        "aPTT_ratio": "aPTT ratio",
        "Xa_inhibition": "Xa inhibition",
    }

    @override
    def models(self) -> dict[str, AbstractModel | Path]:
        return {
            "model": AbstractModel(
                source=MODEL_BASE_PATH / "rivaroxaban_coagulation.xml",
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    def default_changes(self) -> dict[str, Quantity | float]:
        """Changes of every simulation of the experiments, none at the moment.

        The parameters are part of the model and are not changed in the
        simulations: a parameter which is fitted must not be set by a simulation,
        it would override the value of the fit.
        """
        return {}

    @override
    def tasks(self) -> dict[str, Task]:
        return {
            f"task_{key}": Task(model="model", simulation=key)
            for key in self.simulations()
        }

    @override
    def data(self) -> dict[str, Data]:
        self.add_selections_data(
            selections=[
                "time",
                "[Cve_riv]",
                "PT",
                "PT_change",
                "PT_ratio",
                "aPTT",
                "aPTT_change",
                "aPTT_ratio",
                "Xa_inhibition",
            ]
        )
        return {}
