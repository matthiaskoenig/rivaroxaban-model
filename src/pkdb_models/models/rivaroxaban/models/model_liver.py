"""Liver model for rivaroxaban."""

import numpy as np
from sbmlutils.factory import (
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    Reaction,
    Species,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.rivaroxaban.models import annotations, templates
from pkdb_models.models.rivaroxaban.models.fits import FITS


class U(templates.U):
    """Units of the liver model."""


mid = "rivaroxaban_liver"

_m = Model(
    sid=mid,
    name="Model for hepatic rivaroxaban metabolism.",
    notes="""
    Model for hepatic rivaroxaban metabolism.

    Rivaroxaban is metabolized in the liver into a variety of metabolites,
    which are lumped (`rx`): oxidation by CYP3A4 (`f_cyp3a4`) and CYP2J2 and
    CYP independent hydrolysis (Weinz2009, Mueck2014). The metabolites leave
    the liver into the plasma, from which they are excreted in the urine, and
    with the bile into the feces.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),  # Liver
        (BQB.OCCURS_IN, "bto/BTO:0000759"),  # Liver
        (BQB.OCCURS_IN, "ncit/C12392"),  # Liver
        (BQB.HAS_PROPERTY, "ncit/C79371"),  # Pharmacokinetics: Metabolism
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
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True,
    ),
    Compartment(
        "Vmem",
        value=np.nan,
        unit=U.m2,
        name="membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["basolateral"],
        spatialDimensions=2,
        port=True,
    ),
    Compartment(
        "Vbi",
        1.0,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
    ),
    Compartment(
        "Vlumen",
        1.0,
        name="gut lumen",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu"],
        port=True,
    ),
]

_m.species = [
    Species(
        "riv_ext",
        name="rivaroxaban (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "riv",
        name="rivaroxaban (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
    ),
    Species(
        "rx_ext",
        name="rivaroxaban metabolites (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
    Species(
        "rx",
        name="rivaroxaban metabolites (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
    ),
    Species(
        "rx_bi",
        initialConcentration=0.0,
        name="rivaroxaban metabolites (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
    ),
    Species(
        "rx_lumen",
        initialConcentration=0.0,
        name="rivaroxaban metabolites (gut)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "RXEX_BI_k",
        **FITS.parameter("LI__RXEX_BI_k", unit=U.per_min),
        name="rate for rivaroxaban metabolite export in bile",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        Biliary excretion of the metabolites: of an oral solution of 10 mg
        [14C]rivaroxaban 65.7 % of the dose are recovered in the urine (36 %
        unchanged) and 28.0 % in the feces (7 % unchanged), i.e. about 30 % of
        the dose are metabolites in the urine and 21 % in the feces
        (Weinz2009).
        """,
    ),
    Parameter(
        "f_cyp3a4",
        1.0,
        unit=U.dimensionless,
        name="CYP3A4 activity",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Activity of CYP3A4 relative to normal: 1.0 normal, <1.0 inhibition,
        >1.0 induction; set by the experiments of the coadministrations with
        documented factors, not fitted.
        """,
    ),
    Parameter(
        "fm_cyp3a4",
        0.4,
        unit=U.dimensionless,
        name="fraction of the hepatic metabolism by CYP3A4",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""
        Rivaroxaban is metabolized by CYP3A4 (about 18 % of the dose), CYP2J2
        (about 14 %) and CYP independent hydrolysis (about 14 %) (Mueck2014),
        so CYP3A4 is about 18/46 = 0.4 of the hepatic metabolism.
        """,
    ),
]
_m.reactions = [
    Reaction(
        sid="RIVIM",
        name="rivaroxaban import (RIVIM)",
        equation="riv_ext <-> riv",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RIVIM_k",
                1000.0,
                U.per_min,
                name="rate of rivaroxaban import",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("RIVIM_k * Vli * (riv_ext - riv)", U.mmole_per_min),
        notes="""Assuming fast equilibration.""",
    ),
    Reaction(
        sid="RIV2RX",
        name="rivaroxaban conversion (RIV2RX)",
        equation="riv -> rx",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "RIV2RX_k",
                **FITS.parameter("LI__RIV2RX_k", unit=U.per_min),
                name="rate of rivaroxaban conversion",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                Intrinsic hepatic clearance per liver volume; the total
                clearance of rivaroxaban is about 10 l/hr, of which about 6-7
                l/hr are non-renal (Mueck2014).
                """,
            ),
        ],
        formula=(
            "(1 dimensionless - fm_cyp3a4 + fm_cyp3a4 * f_cyp3a4)"
            " * RIV2RX_k * Vli * riv",
            U.mmole_per_min,
        ),
        notes="""
        Oxidative degradation of the morpholinone moiety by CYP3A4 and CYP2J2
        and CYP independent hydrolysis of the amide bonds into a variety of
        metabolites, which are lumped (Weinz2009); the part of CYP3A4 is
        changed by `f_cyp3a4`.
        """,
    ),
    Reaction(
        sid="RXEX",
        name="rivaroxaban metabolite export (RXEX)",
        equation="rx <-> rx_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEX_k",
                1000.0,
                U.per_min,
                name="rate of rivaroxaban metabolite export",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("RXEX_k * Vli * (rx - rx_ext)", U.mmole_per_min),
        notes="""Assuming fast equilibration.""",
    ),
    Reaction(
        "RXEX_BI",
        name="rivaroxaban metabolite bile export",
        equation="rx -> rx_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vmem",
        formula=(
            "RXEX_BI_k * Vli * rx",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "RXEXEHC",
        name="rivaroxaban metabolite enterohepatic circulation",
        equation="rx_bi -> rx_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        formula=(
            "RXEX_BI",
            U.mmole_per_min,
        ),
    ),
]

model_liver = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
