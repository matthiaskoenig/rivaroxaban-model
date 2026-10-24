"""Model for coagulation factors and readouts."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    Species,
    UnitDefinition,
    ValidationOptions,
    create_model,
)
from sbmlutils.metadata import (
    SBO,
)

from pkdb_models.models.rivaroxaban.models import annotations, templates
from pkdb_models.models.rivaroxaban.models.fits import FITS


class U(templates.U):
    """Units of the coagulation model."""

    per_mM = UnitDefinition("per_mM", "liter/mmole")


mid = "rivaroxaban_coagulation"

_m = Model(
    sid=mid,
    name="Model for coagulation factors and readouts.",
    notes="""
    Model for coagulation factors and readouts.

    The effects depend on the total concentration of rivaroxaban in the venous
    plasma `Cve_riv`, as the assays measure the total drug in a plasma sample.
    The prothrombin time (PT) is linear in the concentration (Neoplastin, up
    to the concentrations of 80 mg, Kubitza2005a, Mueck2014); the activated
    partial thromboplastin time (aPTT) and the inhibition of factor Xa follow
    an Emax model of the concentration, the PT and the aPTT relative to their
    baseline (`PT_ref`, `aPTT_ref`), which an experiment sets to the baseline
    of its subjects.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model,
)

_m.compartments = [
    Compartment(
        "Vplasma",
        value=5.0,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True,
    ),
]

_m.species = [
    Species(
        "Cve_riv",
        name="rivaroxaban",
        initialConcentration=0,
        compartment="Vplasma",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
    )
]

names = {
    "PT": "prothrombin time",
    "aPTT": "activated partial thromboplastin time",
    "Xa_inhibition": "inhibition of factor Xa",
}
reference_values = {
    "PT": 12.5,  # [s]
    "aPTT": 28.4,  # [s]
    "Xa_inhibition": 0,  # [-] not inhibited
}
unit_values = {
    "PT": U.second,
    "aPTT": U.second,
    "Xa_inhibition": U.dimensionless,
}

# Xa, PT aPTT
for sid in ["PT", "aPTT"]:
    _m.parameters.extend(
        [
            Parameter(
                sid=f"{sid}_ref",
                name=f"{names[sid]} reference",
                value=reference_values[sid],
                unit=unit_values[sid],
                constant=True,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        ]
    )

for sid in ["PT", "aPTT", "Xa_inhibition"]:
    _m.parameters.extend(
        [
            Parameter(
                sid=sid,
                name=names[sid],
                value=reference_values[sid],
                unit=unit_values[sid],
                constant=False,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            )
        ]
    )

for sid in ["PT", "aPTT"]:
    # absolute change
    _m.parameters.append(
        Parameter(
            f"{sid}_change",
            name=f"{names[sid]} change",
            value=np.nan,
            unit=U.second,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"Absolute change to baseline {names[sid]}",
        )
    )
    _m.rules.append(AssignmentRule(f"{sid}_change", f"{sid}-{sid}_ref", unit=U.second))
    # ratio to baseline
    _m.parameters.append(
        Parameter(
            f"{sid}_ratio",
            name=f"{names[sid]} ratio",
            value=np.nan,
            unit=U.dimensionless,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"Ratio relative to baseline {names[sid]}",
        )
    )
    _m.rules.append(
        AssignmentRule(f"{sid}_ratio", f"{sid}/{sid}_ref", unit=U.dimensionless)
    )


# Effect of rivaroxaban on PT, aPTT and Xa_inhibition
_m.parameters.extend(
    [
        Parameter(
            sid="PT_slope",
            name="relative increase of the PT per rivaroxaban concentration",
            **FITS.parameter("PT_slope", unit=U.per_mM),
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The PT is linear in the total concentration of rivaroxaban in
            plasma; the PT of 10 mg in the fasted state rises to about 1.44
            times its baseline at a Cmax of about 0.26 µM (Neoplastin Plus,
            Kubitza2006), the less sensitive reagents of the other studies give
            smaller increases.
            """,
        ),
        Parameter(
            sid="Emax_aPTT",
            name="maximum effect aPTT",
            **FITS.parameter("Emax_aPTT", unit=U.dimensionless),
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            sid="EC50_riv_aPTT",
            name="EC50 rivaroxaban aPTT",
            **FITS.parameter("EC50_riv_aPTT", unit=U.mM),
            sboTerm=SBO.EC50,
            notes="Total concentration of rivaroxaban in plasma.",
        ),
        Parameter(
            sid="Emax_Xa",
            name="maximum effect Xa inhibition",
            **FITS.parameter("Emax_Xa", unit=U.dimensionless),
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            sid="EC50_riv_Xa",
            name="EC50 rivaroxaban Xa inhibition",
            **FITS.parameter("EC50_riv_Xa", unit=U.mM),
            sboTerm=SBO.EC50,
            notes="Total concentration of rivaroxaban in plasma.",
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "PT",
            "PT_ref * (1 dimensionless + PT_slope * Cve_riv)",
            unit=U.second,
        ),
        AssignmentRule(
            "aPTT",
            "aPTT_ref * (1 dimensionless + Emax_aPTT * Cve_riv / (Cve_riv + EC50_riv_aPTT))",
            unit=U.second,
        ),
        AssignmentRule(
            "Xa_inhibition",
            "Emax_Xa * Cve_riv / (Cve_riv + EC50_riv_Xa)",
            unit=U.dimensionless,
        ),
    ]
)

model_coagulation = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_coagulation,
        filepath=MODEL_BASE_PATH / f"{model_coagulation.sid}.xml",
        sbml_level=3,
        sbml_version=2,
        validation_options=ValidationOptions(units_consistency=True),
    )
    # create differential equations
    md_path = MODEL_BASE_PATH / f"{model_coagulation.sid}.md"
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(md_path)
