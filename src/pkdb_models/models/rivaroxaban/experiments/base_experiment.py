"""
Reusable functionality for multiple simulation experiments.
"""
import pandas as pd
from collections import namedtuple
from typing import Dict
from pkdb_models.models.rivaroxaban import MODEL_PATH
from pkdb_models.models.rivaroxaban.rivaroxaban_pk import (
    calculate_rivaroxaban_pk,
    calculate_rivaroxaban_pd,
)
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "riv rx")


class RivaroxabanSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
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

    labels: Dict[str, str] = {
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

        "PT": "PT",  # "prothrombin time",
        "PT_change": "PT change", # "prothrombin time change",
        "PT_ratio": "PT ratio",  # "prothrombin time ratio",
        "aPTT": "aPTT",
        "aPTT_change": "aPTT change",
        "aPTT_ratio": "aPTT ratio",
        "Xa_inhibition": "Xa inhibition",
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

    units: Dict[str, str] = {
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

        "PT": "s",
        "PT_change": "s",
        "PT_ratio": "dimensionless",
        "aPTT": "s",
        "aPTT_change": "s",
        "aPTT_ratio": "dimensionless",
        "Xa_inhibition": "dimensionless",
    }

    # ----------- Fasting/food -----
    # food changes the fraction absorbed
    fasting_map = {  # GU__F_riv_abs
        "not reported": 0.82,  # assuming fasted state if nothing is reported
        "fasted": 0.82,
        "fed": 1.0,
    }
    fasting_colors = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 50.0 / 101.0,  # 0.5
        "Moderate renal impairment": 35.0 / 101.0,  # 0.35
        "Severe renal impairment": 20.0 / 101.0,  # 0.20
        # "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        # "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",  # CPT C
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {

            # Optimization results: 20250717_161243__4822f/RIVAROXABAN_LSQ_PK
            # >>> !Optimal parameter 'LI__RXEX_BI_k' within 5% of upper bound! <<<
            # 'GU__Ka_dis_riv': Q_(0.15960999342060214, '1/hr'),  # [0.001 - 100]
            # 'GU__RIVABS_k': Q_(0.052027850297051675, '1/min'),  # [0.001 - 10]
            # 'LI__RIV2RX_Vmax': Q_(0.06401213063939644, '1/min'),  # [0.001 - 100]
            # 'LI__RXEX_BI_k': Q_(9.99999999999998e-05, '1/min'),  # [1e-06 - 0.0001]
            # 'KI__RIVEX_k': Q_(0.16139673072970187, '1/min'),  # [0.0001 - 10]


            # pharmacodynamic
            # Optimization results: 20250717_162628__71845/RIVAROXABAN_LSQ_PD
            # 'Emax_PT': Q_(3.5530919743774545, 'dimensionless'),  # [0.1 - 10]
            # 'EC50_riv_PT': Q_(0.0035833544566150222, 'mM'),  # [1e-07 - 0.01]
            # 'Emax_aPTT': Q_(0.9500194161813726, 'dimensionless'),  # [0.1 - 10]
            # 'EC50_riv_aPTT': Q_(0.00040891827752344355, 'mM'),  # [1e-07 - 0.01]
            # 'Emax_Xa': Q_(0.6860934265985122, 'dimensionless'),  # [0.5 - 1.0]
            # 'EC50_riv_Xa': Q_(0.000295895101103941, 'mM'),  # [1e-07 - 0.01]
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return RivaroxabanSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
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

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',

                "PT",
                "PT_change",
                "PT_ratio",
                "aPTT",
                "aPTT_change",
                "aPTT_ratio",
                "Xa_inhibition",

                "PODOSE_riv",
                "GU__F_riv_abs",
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            riv=self.Q_(435.88, "g/mole"),
            rx=self.Q_(435.88, "g/mole"),  # FIXME
        )

    # --- Pharmacokinetic parameters ---
    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        "cl_renal": "Renal clearance",
        "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "µmole/l*hr",
        "cl": "ml/min",
        "cl_renal": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
    }

    def calculate_rivaroxaban_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pk parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rivaroxaban_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rivaroxaban_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs

    def calculate_rivaroxaban_pd(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pd parameters for simulations (scans)"""
       pd_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rivaroxaban_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_rivaroxaban_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       return pd_dfs