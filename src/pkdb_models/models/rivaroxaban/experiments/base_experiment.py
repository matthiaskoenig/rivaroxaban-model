"""Base class of the rivaroxaban simulation experiments."""

from collections.abc import Mapping, Sequence
from typing import ClassVar, NamedTuple

import pandas as pd
from pkpdsim.experiment import PKPDSimulationExperiment
from pkpdsim.impairment import RENAL_FUNCTION_GFR
from sbmlsim import Q
from sbmlsim.units import Quantity

from pkdb_models.models.rivaroxaban import MODEL
from pkdb_models.models.rivaroxaban.rivaroxaban_pk import (
    calculate_rivaroxaban_pd,
    calculate_rivaroxaban_pk,
)


class MolecularWeights(NamedTuple):
    """Molecular weights for the conversion of the data."""

    riv: Quantity
    rx: Quantity


class RivaroxabanSimulationExperiment(PKPDSimulationExperiment):
    """Base class of the rivaroxaban simulation experiments."""

    pkpd_model = MODEL
    selections = (
        "time",
        "[Cve_riv]",
        "[Cve_rx]",
        "[Cve_riv_rx]",
        "Aurine_riv",
        "Aurine_rx",
        "Aurine_riv_rx",
        "Afeces_riv",
        "Afeces_rx",
        "Afeces_riv_rx",
        "Afeces_urine_riv_rx",
        "KI__RIVEX",
        "KI__GFR",
        # cases
        "KI__f_renal_function",
        "f_cirrhosis",
        "LI__f_cyp3a4",
        "KI__f_pgp",
        "GU__fed",
        "PT",
        "PT_change",
        "PT_ratio",
        "aPTT",
        "aPTT_change",
        "aPTT_ratio",
        "Xa_inhibition",
        "PODOSE_riv",
        "POSOLUTION_riv",
        "GU__POSTOMACH_riv",
        "GU__POINTESTINE_riv",
        "IVDOSE_riv",
    )

    font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 22}
    scan_font: ClassVar[dict[str, str | int]] = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "time"
    label_riv = "rivaroxaban"
    label_rx = "rivaroxaban metabolites"
    label_riv_rx = "RIV + RX"

    label_riv_plasma = label_riv + " plasma"
    label_rx_plasma = label_rx + "\nplasma"
    label_riv_rx_plasma = label_riv_rx + "\nplasma"
    label_riv_urine = label_riv + " urine"
    label_rx_urine = label_rx + "\nurine"
    label_riv_rx_urine = label_riv_rx + "\nurine"
    label_riv_feces = label_riv + "\nfeces"
    label_rx_feces = label_rx + "\nfeces"
    label_riv_rx_feces = label_riv_rx + "\nfeces"

    labels: ClassVar[Mapping[str, str]] = {
        "time": "time",
        "[Cve_riv]": label_riv_plasma,
        "[Cve_rx]": label_rx_plasma,
        "[Cve_riv_rx]": label_riv_rx_plasma,
        "Aurine_riv": label_riv_urine,
        "Aurine_rx": label_rx_urine,
        "Aurine_riv_rx": label_riv_rx_urine,
        "Afeces_riv": label_riv_feces,
        "Afeces_rx": label_rx_feces,
        "Afeces_riv_rx": label_riv_rx_feces,
        "KI__RIVEX": label_riv + " excretion\nurine",
        "KI__GFR": "GFR",
        "PT": "PT",  # "prothrombin time",
        "PT_change": "PT change",  # "prothrombin time change",
        "PT_ratio": "PT ratio",  # "prothrombin time ratio",
        "aPTT": "aPTT",
        "aPTT_change": "aPTT change",
        "aPTT_ratio": "aPTT ratio",
        "Xa_inhibition": "Xa inhibition",
        "PODOSE_riv": "oral dose\n(not applied)",
        "POSOLUTION_riv": "oral solution\n(not applied)",
        "GU__POSTOMACH_riv": "oral dose\n(stomach)",
        "GU__POINTESTINE_riv": "oral dose\n(not dissolved)",
        "IVDOSE_riv": "intravenous dose\n(not injected)",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "µM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_riv = unit_metabolite
    unit_rx = unit_metabolite
    unit_riv_rx = unit_metabolite

    unit_riv_urine = unit_metabolite_urine
    unit_rx_urine = unit_metabolite_urine
    unit_riv_rx_urine = unit_metabolite_urine
    unit_riv_feces = unit_metabolite_feces
    unit_rx_feces = unit_metabolite_feces
    unit_riv_rx_feces = unit_metabolite_feces

    units: ClassVar[Mapping[str, str]] = {
        "time": unit_time,
        "[Cve_riv]": unit_riv,
        "[Cve_rx]": unit_rx,
        "[Cve_riv_rx]": unit_riv_rx,
        "Aurine_riv": unit_riv_urine,
        "Aurine_rx": unit_rx_urine,
        "Aurine_riv_rx": unit_riv_rx_urine,
        "Afeces_riv": unit_riv_feces,
        "Afeces_rx": unit_rx_feces,
        "Afeces_riv_rx": unit_riv_rx_feces,
        "Afeces_urine_riv_rx": unit_riv_rx_feces,
        "KI__RIVEX": "µmole/hr",
        "KI__GFR": "ml/min",
        "PT": "s",
        "PT_change": "s",
        "PT_ratio": "dimensionless",
        "aPTT": "s",
        "aPTT_change": "s",
        "aPTT_ratio": "dimensionless",
        "Xa_inhibition": "dimensionless",
        "PODOSE_riv": "mg",
        "POSOLUTION_riv": "mg",
        "GU__POSTOMACH_riv": "mg",
        "GU__POINTESTINE_riv": "mg",
        "IVDOSE_riv": "mg",
    }

    # ----------- Fasting/food -----
    #: the fed state of the intestine (`GU__fed`): food increases the
    #: solubility and slows the gastric emptying, see the intestine model
    fasting_map: ClassVar[dict[str, float]] = {
        "not reported": 0.0,  # assuming fasted state if nothing is reported
        "fasted": 0.0,
        "fed": 1.0,
    }
    fasting_colors: ClassVar[dict[str, str]] = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal map --------------
    # the classes of RENAL_FUNCTION_GFR without the end stage renal disease;
    # the renal function is relative to the GFR of 101 ml/min, a measured
    # creatinine clearance of X ml/min is the renal function X/101
    renal_map: ClassVar[Mapping[str, float]] = {
        key: value
        for key, value in RENAL_FUNCTION_GFR.items()
        if key != "End stage renal disease"
    }

    #: molecular weights of rivaroxaban and its metabolites; the metabolites
    #: are lumped in rivaroxaban equivalents, as the radioactivity of
    #: [14C]rivaroxaban (Weinz2009), the only data of the metabolites
    Mr: ClassVar[MolecularWeights] = MolecularWeights(
        riv=Q(435.88, "g/mole"),
        rx=Q(435.88, "g/mole"),
    )

    def calculate_rivaroxaban_pk(
        self, scans: Sequence[str] | None = None, across: str | None = None
    ) -> dict[str, pd.DataFrame]:
        """Calculate the pharmacokinetic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.
            across: the dimension of the draws of the fitted parameters,
                `None` for scans without draws.

        Returns:
            The parameters of every scan by simulation key.
        """
        keys = scans if scans else self._simulations
        return {
            key: calculate_rivaroxaban_pk(
                experiment=self, res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }

    def calculate_rivaroxaban_pd(
        self, scans: Sequence[str] | None = None, across: str | None = None
    ) -> dict[str, dict[str, pd.DataFrame]]:
        """Calculate the pharmacodynamic parameters of the scans.

        Args:
            scans: keys of the simulations, all simulations by default.
            across: the dimension of the draws of the fitted parameters,
                `None` for scans without draws.

        Returns:
            The parameters of every scan by simulation key and effect.
        """
        keys = scans if scans else self._simulations
        return {
            key: calculate_rivaroxaban_pd(
                res=self.results[f"task_{key}"], across=across
            )
            for key in keys
        }
