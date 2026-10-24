"""Kidney model for rivaroxaban."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    Reaction,
    Species,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.rivaroxaban.models import annotations, templates
from pkdb_models.models.rivaroxaban.models.fits import FITS


class U(templates.U):
    """Units of the kidney model."""

    l_per_ml = UnitDefinition("l_per_ml", "l/ml")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")


mid = "rivaroxaban_kidney"

_m = Model(
    sid=mid,
    name="Model for renal rivaroxaban excretion.",
    notes="""
    Model for renal rivaroxaban excretion.

    The renal clearance of rivaroxaban is the filtration of its unbound
    fraction, `GFR * fu_riv`, times `RIVEX_f` for the active secretion by P-gp
    and BCRP (`f_pgp`); the glomerular filtration rate is the renal function
    times the healthy GFR of 101 ml/min. About 36 % of an oral dose of 10 mg is
    excreted unchanged in the urine (Weinz2009, Kubitza2013).
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),
        (BQB.OCCURS_IN, "ncit/C12415"),
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True,
    ),
    Compartment(
        "Vki",
        value=0.3,
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True,
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),
]

_m.species = [
    Species(
        "riv_ext",
        name="rivaroxaban (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "riv_urine",
        name="rivaroxaban (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "rx_ext",
        name="rivaroxaban metabolites (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
    Species(
        "rx_urine",
        name="rivaroxaban metabolites (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
]

# glomerular filtration rate and kidney function
_m.parameters.extend(
    [
        Parameter(
            "f_renal_function",
            name="renal function",
            value=1.0,
            unit=U.dimensionless,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            >1.0: increased kidney function
            1.0: normal kidney function (healthy control)
            <1.0: decreased kidney function

            Renal function relative to the healthy glomerular filtration rate
            GFR_base, which is 101 ml/min at 75 kg; a measured creatinine
            clearance of X ml/min is the renal function X/101.
            """,
        ),
        Parameter(
            "GFR_base",
            101,
            U.ml_per_min,
            name="glomerular filtration rate (healthy)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
        ),
        Parameter(
            "GFR",
            np.nan,
            U.ml_per_min,
            constant=False,
            name="glomerular filtration rate",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "fu_riv",
            0.079,
            U.dimensionless,
            constant=True,
            name="fraction unbound in plasma rivaroxaban",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            7.9 % unbound in the plasma of healthy subjects (Kubitza2013), 92-95 %
            bound to plasma proteins, mainly albumin.
            """,
            port=True,
        ),
        Parameter(
            "RIVEX_f",
            **FITS.parameter("KI__RIVEX_f", unit=U.dimensionless),
            name="renal clearance of unbound rivaroxaban relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of rivaroxaban is RIVEX_f times the filtration
            of its unbound fraction, fu_riv * GFR: 1 is filtration alone, above
            1 a net secretion. Rivaroxaban is secreted by P-gp and BCRP, about
            two thirds of the renal clearance of 2.4-3.1 l/hr are active
            (Mueck2014, Kubitza2013); 36-42 % of a small oral dose are
            excreted unchanged in the urine (Weinz2009, Kubitza2005a), i.e.
            RIVEX_f about 6-8.
            """,
        ),
        Parameter(
            "f_pgp",
            1.0,
            U.dimensionless,
            name="P-gp/BCRP activity",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
            notes="""
            Activity of the renal efflux transporters (P-gp, BCRP) relative to
            normal: 1.0 normal, <1.0 inhibition, >1.0 induction. It scales the
            secretion, `(RIVEX_f - 1)`, not the filtration; set by the
            experiments of the coadministrations with documented factors, not
            fitted.
            """,
        ),
        Parameter(
            "RXEX_f",
            7.0,
            U.dimensionless,
            name="renal clearance of the rivaroxaban metabolites relative to the GFR",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The renal clearance of the metabolites is RXEX_f times the GFR; the
            metabolites are about 11 % of the radioactivity in plasma after
            [14C]rivaroxaban (Weinz2009), so they are cleared faster than they
            are filtered.
            """,
        ),
    ]
)
_m.rules.append(
    AssignmentRule(
        "GFR",
        "f_renal_function * GFR_base",
        unit=U.ml_per_min,
        name="glomerular filtration rate",
    ),
)

_m.reactions = [
    Reaction(
        sid="RIVEX",
        name="rivaroxaban excretion (RIVEX)",
        equation="riv_ext -> riv_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        # [ml/min] * [l/ml] * [-] * [-] * [mmole/l]
        formula=(
            "GFR * 0.001 l_per_ml * fu_riv"
            " * (1 dimensionless + (RIVEX_f - 1 dimensionless) * f_pgp) * riv_ext",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="RXEX",
        name="rivaroxaban metabolites excretion (RXEX)",
        equation="rx_ext -> rx_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "GFR * 0.001 l_per_ml * RXEX_f * rx_ext",
            U.mmole_per_min,
        ),
    ),
]

model_kidney = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
