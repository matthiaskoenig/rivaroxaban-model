"""Sensitivity analysis of the rivaroxaban model.

The outputs are the pharmacokinetic parameters of a single oral dose of 10 mg
(the non-compartmental analysis of `pkpdutils`) of rivaroxaban and its
metabolites and the maxima of the pharmacodynamics. The parameters are the
parameters of the model with bounds of +-15 % around their value, the fitted
parameters with the bounds of the fits. The analyses are run as a script, the
local and the sampling analysis by default:

    python -m pkdb_models.models.rivaroxaban.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.rivaroxaban.sensitivity.sensitivity_analysis local

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which its simulation or one of its groups
sets (`pkpdsim.sensitivity.free_parameters`): every analysis leaves out
`KI__f_renal_function`, which the simulation sets (`f_cirrhosis`, 0 in the
model, is no parameter of the analyses), `sampling`, over the eight groups,
also `GU__fed` (0 in the model, no parameter of the analyses either). The environment variable `PKPDSIM_SENSITIVITY_N` replaces
the numbers of samples of the global analyses for a quick run, see
`pkpdsim.sensitivity.sample_size`.
"""

from collections.abc import Mapping
from typing import Any, NamedTuple

import numpy as np
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from sbmlsim import Q
from sbmlsim.simulation import PK, Custom, Formula, ObservableKind, Simulation

from pkdb_models.models.rivaroxaban import MODEL, MODEL_PATH
from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
    RivaroxabanSimulationExperiment,
)
from pkdb_models.models.rivaroxaban.fitting.parameters import (
    parameters_pd,
    parameters_pk,
)

#: oral dose of the analysis [mg]
DOSE_RIV = 10
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15

#: outputs of the non-compartmental analysis of pkpdutils -> unit of the output
PK_OUTPUTS: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "thalf": "min",
    "vz_f": "l",
    "cl_f": "ml/min",
    "lambda_z": "1/min",
}
#: the outputs of the metabolites, which have no dose of their own
PK_OUTPUTS_RX: list[str] = ["auc_inf_obs", "cmax", "thalf"]
#: limit of quantification relative to the maximal concentration: the
#: concentrations of the 20 days which are smaller are the noise of the
#: integrator, they would be the terminal phase
LLOQ_FRACTION = 1e-8
#: the outputs whose unit has the amount of the dose: the dose is a mass
#: (`PODOSE_riv` in mg), so they are divided by the molecular weight (its value
#: at the start, a value per simulation like the parameter of the analysis)
PK_PER_DOSE_MASS: frozenset[str] = frozenset({"vz_f", "cl_f"})
#: the labels of the pharmacokinetic parameters in the names of the outputs
PK_LABELS: dict[str, str] = {
    "auc_inf_obs": "AUC∞",
    "cmax": "Cmax",
    "thalf": "Half-life",
    "vz_f": "Vz/F",
    "cl_f": "CL/F",
    "lambda_z": "λz",
}
#: the pharmacodynamic outputs, their maximum is the output -> name, unit
PD_OUTPUTS: dict[str, tuple[str, str]] = {
    "PT": ("max PT", "s"),
    "aPTT": ("max aPTT", "s"),
    "Xa_inhibition": ("max Xa inhibition", "dimensionless"),
}


def quantified(time: np.ndarray, values: dict[str, Any]) -> np.ndarray:
    """Get a concentration above its limit of quantification, missing below it.

    The limit is `LLOQ_FRACTION` of the maximum of the timecourse, as the
    limit of quantification of the old analysis (`Timecourse.lloq`), whose
    values below it the non-compartmental analysis leaves out.

    Args:
        time: the time points of the simulation.
        values: the concentration, the one symbol of the observable.

    Returns:
        The concentration, `NaN` below the limit.
    """
    (concentration,) = values.values()
    c = np.asarray(concentration, dtype=float)
    return np.where(c < LLOQ_FRACTION * np.max(c), np.nan, c)


class Substance(NamedTuple):
    """A substance of the pharmacokinetic outputs, its concentration in plasma."""

    #: the concentration, which the non-compartmental analysis reads
    selection: str
    #: the name of the substance in the names of the outputs
    name: str
    #: the id of the non-compartmental analysis, no id of the model
    pk_id: str
    #: the outputs of the analysis -> unit
    outputs: dict[str, str]
    #: the dose, of rivaroxaban only
    dose: str | None

    @property
    def quantified_id(self) -> str:
        """Get the id of the concentration above the limit of quantification."""
        return f"{self.selection.strip('[]')}_quantified"


#: the substances of the pharmacokinetic outputs
SUBSTANCES: list[Substance] = [
    Substance("[Cve_riv]", "RIV", "pk_riv", PK_OUTPUTS, "PODOSE_riv"),
    Substance(
        "[Cve_rx]",
        "RX",
        "pk_rx",
        {key: PK_OUTPUTS[key] for key in PK_OUTPUTS_RX},
        None,
    ),
]

#: the outputs: the id of a pharmacokinetic output names the concentration with
#: its brackets (`[Cve_riv]_cmax`), as the tables of sbmlsim, its
#: observable has the id without them (`Cve_riv_cmax`)
OUTPUTS: list[SensitivityOutput] = [
    *(
        SensitivityOutput(
            f"{s.selection}_{key}",
            f"{s.name} {PK_LABELS[key]}",
            f"{s.selection.strip('[]')}_{key}",
        )
        for s in SUBSTANCES
        for key in s.outputs
    ),
    *(SensitivityOutput(f"{sid}_max", name) for sid, (name, _) in PD_OUTPUTS.items()),
]


def _pk_formula(substance: Substance, key: str) -> str:
    """Get the formula of a pharmacokinetic output of a substance."""
    if key in PK_PER_DOSE_MASS:
        return f"{substance.pk_id}.{key} / at(Mr_riv, 0)"
    return f"{substance.pk_id}.{key}"


#: the observables of a simulation, the outputs with their units
OBSERVABLES = [
    *(
        Custom(
            s.quantified_id,
            quantified,
            "mM",
            symbols=[s.selection],
            kind=ObservableKind.TIMECOURSE,
        )
        for s in SUBSTANCES
    ),
    *(
        PK(
            s.pk_id,
            s.quantified_id,
            dose=s.dose,
            route="oral" if s.dose else None,
            parameters=list(s.outputs),
        )
        for s in SUBSTANCES
    ),
    *(
        Formula(f"{s.selection.strip('[]')}_{key}", _pk_formula(s, key), unit=unit)
        for s in SUBSTANCES
        for key, unit in s.outputs.items()
    ),
    *(
        Formula(f"{sid}_max", f"max({sid})", unit=unit)
        for sid, (_, unit) in PD_OUTPUTS.items()
    ),
]


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: a single oral dose of rivaroxaban.

    The simulation sets the healthy liver and kidney, which a group changes.
    """
    return Simulation(
        end=20 * 24 * 60,
        time_unit="min",
        preinit_changes={
            "PODOSE_riv": Q(DOSE_RIV, "mg"),
            "f_cirrhosis": Q(0, "dimensionless"),
            "KI__f_renal_function": Q(1.0, "dimensionless"),
            **changes,
        },
    )


#: the groups the analyses are run for: the control, the stages of renal
#: impairment and of cirrhosis of the simulation experiments (`renal_map`,
#: `cirrhosis_map`) and the fed state
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {
                    "KI__f_renal_function": RivaroxabanSimulationExperiment.renal_map[
                        name
                    ]
                }
            ),
            RivaroxabanSimulationExperiment.renal_colors[name],
        )
        for uid, name in [
            ("mildRI", "Mild renal impairment"),
            ("modRI", "Moderate renal impairment"),
            ("sevRI", "Severe renal impairment"),
        ]
    ),
    *(
        SensitivityGroup(
            uid,
            f"{name} ({uid})",
            sensitivity_simulation(
                {"f_cirrhosis": RivaroxabanSimulationExperiment.cirrhosis_map[name]}
            ),
            RivaroxabanSimulationExperiment.cirrhosis_colors[name],
        )
        for uid, name in [
            ("CPT A", "Mild cirrhosis"),
            ("CPT B", "Moderate cirrhosis"),
            ("CPT C", "Severe cirrhosis"),
        ]
    ),
    SensitivityGroup(
        "fed", "Fed state", sensitivity_simulation({"GU__fed": 1}), "tab:red"
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds and references."""
    parameters = parameters_from_sbml(
        MODEL_PATH,
        exclude_ids={
            # conversion factors
            "conversion_min_per_day",
            # molecular weights
            "Mr_riv",
            "Mr_rx",
            # unchangeable values
            "FQlu",
            "FVhv",
            "FVpo",
            # fast transports
            "LI__RIVIM_k",
            "LI__RXEX_k",
            "GU__Ka_application_riv",
            # dosing parameters
            "PODOSE_riv",
            "POSOLUTION_riv",
            "IVDOSE_riv",
            "ti_riv",
            "Ki_riv",
            "Ri_riv",
            # unused volumes
            "Vurine",
            "Vfeces",
            "Vplasma",
            "LI__Vbi",
            # reference values
            "PT_ref",
            "aPTT_ref",
            "BW_ref",
            "KI__GFR_base",
        },
        exclude_na=True,
        exclude_zero=True,
    )

    # the fitted parameters vary within the bounds of the fits
    set_bounds(
        parameters,
        [
            (fp.pid, fp.lower_bound, fp.upper_bound, ParameterType.FIT)
            for fp in [*parameters_pk, *parameters_pd]
        ],
    )
    # the renal function and the fraction absorbed, which cannot exceed 1
    set_bounds(
        parameters,
        [
            (
                "KI__f_renal_function",
                1 - BOUNDS_FRACTION,
                1 + BOUNDS_FRACTION,
                ParameterType.SCALING,
            ),
            ("GU__F_riv_abs", 0.85, 1.0, ParameterType.DATA),
        ],
    )

    # the references of the values from the literature
    references = {
        "HCT": r"\cite{Mondal2025, Fiseha2023}",
        "BW": r"\cite{Ogden2004, Jones2013, Thompson2009, Brown1997}",
        "FQgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FQh": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FQki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVar": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVli": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FVlu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVve": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
        "GU__F_riv_abs": r"\cite{Weinz2009}",
        "KI__fu_riv": r"\cite{Kubitza2013}",
        "LI__fm_cyp3a4": r"\cite{Mueck2014}",
        "KI__f_renal_function": r"\cite{Stevens2024}",
    }
    by_uid = {p.uid: p for p in parameters}
    for uid, reference in references.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS[:1], delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "fast": FastAnalysis(groups=SENSITIVITY_GROUPS[:1], n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local", "sampling"),
)


if __name__ == "__main__":
    SENSITIVITY.main()
