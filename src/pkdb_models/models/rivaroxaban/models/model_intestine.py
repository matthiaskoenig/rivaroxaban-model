"""Rivaroxaban intestine model."""

import numpy as np
from sbmlutils.factory import (
    AssignmentRule,
    Compartment,
    FactoryResult,
    Model,
    Parameter,
    RateRule,
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
    """Units of the intestine model."""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")
    mg_per_hr = UnitDefinition("mg_per_hr", "mg/hr")
    min_per_hr = UnitDefinition("min_per_hr", "min/hr")


_m = Model(
    "rivaroxaban_intestine",
    name="Model for rivaroxaban absorption in the small intestine",
    notes="""
    # Model for rivaroxaban absorption

    An oral dose of a tablet (`PODOSE_riv`) enters the stomach
    (`POSTOMACH_riv`) and is emptied into the intestine in three stages
    (`emptying<k>_riv`, `POINTESTINE_riv`), where it dissolves into the lumen
    (`dissolution_riv`).
    The dissolution of rivaroxaban, a compound of low solubility (BCS class
    II), saturates at a small undissolved amount (`Km_dis_riv`), its capacity
    is `Vmax_dis_riv`; the undissolved drug leaves the absorptive part of the
    intestine with the transit (`ktransit_riv`) and is excreted with the feces
    (`transit_riv`). The fraction of the dose which dissolves therefore falls
    with the dose: the bioavailability is 80-100 % at 10 mg and about 66 % at
    20 mg in the fasted state, the AUC per dose falls from 10 to 80 mg to 40 %
    and the urinary excretion of unchanged drug from 42 % at 1.25 mg to 11 %
    of the dose at 80 mg (Kubitza2005a); the Cmax of 10, 20, 40 and 80 mg
    (141, 173, 234 and 316 µg/l) follows from the capacity of the dissolution.

    Food (`fed`) increases the solubility (`FED_f_solubility`, bile micelles)
    and slows the gastric emptying (`FED_f_emptying`): with food the AUC of 20
    mg increases by 39 % and the Cmax by 76 %, the tmax is delayed by 1.25 hr;
    at 10 mg the AUC increases by 28 % (Kubitza2006).

    An oral solution (`POSOLUTION_riv`) enters the lumen of the intestine
    without the dissolution; the passage of the solution from the stomach is
    lumped into the absorption (tmax 0.5 hr, Kubitza2005a, Weinz2009).

    The dissolved rivaroxaban is absorbed (`RIVABS`); a fraction `1 -
    F_riv_abs` is not absorbed and excreted with the feces (`RIVEXC`): 7 % of
    an oral solution of [14C]rivaroxaban is excreted unchanged with the feces
    (Weinz2009). The metabolites which are excreted with the bile are excreted
    with the feces (`RXEXC`).
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "ncit/C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine
        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vlumen",
        value=1.0,
        name="lumen",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1.0,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["apical"],
    ),
]

_m.species = [
    Species(
        "riv_ext",
        initialConcentration=0.0,
        name="rivaroxaban (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "riv_lumen",
        initialConcentration=0.0,
        name="rivaroxaban (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "riv_feces",
        initialConcentration=0.0,
        name="rivaroxaban (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        port=True,
    ),
    Species(
        "rx_lumen",
        initialConcentration=0.0,
        name="rivaroxaban metabolites (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
    Species(
        "rx_feces",
        initialConcentration=0.0,
        name="rivaroxaban metabolites (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
]

# -------------------------------------
# Absorption
# -------------------------------------
_m.parameters = [
    Parameter(
        "F_riv_abs",
        0.93,
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name="fraction absorbed of the dissolved rivaroxaban",
        notes="""
        Fraction of the dissolved rivaroxaban in the lumen of the intestine
        which is absorbed, the rest `(1 - F_riv_abs)` is excreted with the
        feces. 7 % of an oral solution of 10 mg [14C]rivaroxaban is excreted
        unchanged with the feces (Weinz2009), so 0.93. The dose dependent and
        the food dependent bioavailability of the tablets are the dissolution,
        see `Vmax_dis_riv` and `fed`.
        """,
    ),
    Parameter(
        "RIVABS_k",
        **FITS.parameter("GU__RIVABS_k", unit=U.per_min),
        name="rate of rivaroxaban absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        First order absorption of the dissolved rivaroxaban, including the
        passage of an oral solution from the stomach: the tmax of an oral
        solution is about 0.5 hr (Kubitza2005a, Weinz2009).
        """,
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption_riv",
        value="RIVABS_k * Vlumen * riv_lumen",
        unit=U.mmole_per_min,
        name="absorption rivaroxaban",
    ),
)

_m.reactions = [
    Reaction(
        "RIVABS",
        name="absorption rivaroxaban",
        equation="riv_lumen -> riv_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_riv_abs * absorption_riv", U.mmole_per_min),
    ),
    Reaction(
        sid="RIVEXC",
        name="excretion rivaroxaban (feces)",
        compartment="Vlumen",
        equation="riv_lumen -> riv_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            "(1 dimensionless - F_riv_abs) * absorption_riv",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="RXEXC",
        name="excretion rivaroxaban metabolites (feces)",
        compartment="Vlumen",
        equation="rx_lumen -> rx_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEXC_k",
                1 / (36 * 60),
                unit=U.per_min,
                name="rate of rivaroxaban metabolites excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                Transit of the metabolites which are excreted with the bile
                into the feces, a mean transit time of about 36 hr: the fecal
                radioactivity of [14C]rivaroxaban is recovered within 2-4 days
                (Weinz2009). The rate does not change the amount which is
                excreted.
                """,
            ),
        ],
        formula=(
            "RXEXC_k * Vlumen * rx_lumen",
            U.mmole_per_min,
        ),
    ),
]

# -------------------------------------
# Food
# -------------------------------------
_m.parameters.extend(
    [
        Parameter(
            "fed",
            0.0,
            unit=U.dimensionless,
            name="fed state (0 fasted, 1 fed)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The dose is taken in the fasted state (0) or with a meal (1). The
            effects of food on the solubility and on the gastric emptying are
            `FED_f_solubility` and `FED_f_emptying`.
            """,
        ),
        Parameter(
            "FED_f_solubility",
            **FITS.parameter("GU__FED_f_solubility", unit=U.dimensionless),
            name="fed effect on the solubility",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Factor of the capacity of the dissolution `Vmax_dis_riv` with food
            (bile micelles): food increases the AUC of 20 mg by 39 % and of 10
            mg by 28 % (Kubitza2006), the bioavailability of 20 mg is about 66
            % fasted and about 100 % fed.
            """,
        ),
        Parameter(
            "FED_f_emptying",
            **FITS.parameter("GU__FED_f_emptying", unit=U.dimensionless),
            name="fed effect on the gastric emptying",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Factor of the rate of the gastric emptying `Ka_emptying_riv` with
            food: food delays the tmax by about 1.25 hr (Kubitza2006).
            """,
        ),
        Parameter(
            "f_solubility",
            np.nan,
            unit=U.dimensionless,
            constant=False,
            name="factor of the solubility",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_emptying",
            np.nan,
            unit=U.dimensionless,
            constant=False,
            name="factor of the gastric emptying",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    ]
)
_m.rules.extend(
    [
        AssignmentRule(
            "f_solubility",
            "1 dimensionless + fed * (FED_f_solubility - 1 dimensionless)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "f_emptying",
            "1 dimensionless + fed * (FED_f_emptying - 1 dimensionless)",
            unit=U.dimensionless,
        ),
    ]
)

# -------------------------------------
# Dosing, gastric emptying and dissolution
# -------------------------------------
#: the stages of the gastric emptying: a chain of first order steps, whose
#: delay is distributed narrower than the one of a single step (a lag)
N_EMPTYING = 3

_m.parameters.extend(
    [
        Parameter(
            "PODOSE_riv",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose rivaroxaban [mg]",
            notes="""
            An oral dose of a tablet; it enters the stomach (`POSTOMACH_riv`)
            within minutes, so that a dose of a repeated dosing is added by
            setting `PODOSE_riv` to the dose.
            """,
            port=True,
        ),
        Parameter(
            "POSOLUTION_riv",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose as solution rivaroxaban [mg]",
            notes="""
            An oral dose given dissolved (a solution) enters the lumen of the
            intestine without the dissolution of `PODOSE_riv`; the passage of
            the solution from the stomach is lumped into the absorption.
            """,
            port=True,
        ),
        Parameter(
            "POSTOMACH_riv",
            np.nan,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose in stomach rivaroxaban [mg]",
            notes="""
            The dose in the stomach, the sum of the stages of the gastric
            emptying (`POSTOMACH<k>_riv`).
            """,
        ),
        *(
            Parameter(
                f"POSTOMACH{k}_riv",
                0,
                U.mg,
                constant=False,
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                name=f"oral dose in stomach rivaroxaban, stage {k} [mg]",
            )
            for k in range(1, N_EMPTYING + 1)
        ),
        Parameter(
            "POINTESTINE_riv",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="undissolved oral dose in intestine rivaroxaban [mg]",
            notes="""
            The undissolved rivaroxaban in the absorptive part of the
            intestine.
            """,
        ),
        Parameter(
            "Ka_application_riv",
            1000,
            U.per_hr,
            constant=True,
            sboTerm=SBO.KINETIC_CONSTANT,
            name="Ka [1/hr] application rivaroxaban",
            notes="""Fast application to shift the applied dose into the stomach.""",
        ),
        Parameter(
            "application_riv",
            np.nan,
            U.mmole_per_min,
            constant=False,
            sboTerm=SBO.RATE_OF_REACTION,
            name="PO application rivaroxaban",
            notes="""
            The rate at which the oral dose `PODOSE_riv` enters the stomach
            (`POSTOMACH1_riv`).
            """,
        ),
        Parameter(
            "Ka_emptying_riv",
            **FITS.parameter("GU__Ka_emptying_riv", unit=U.per_hr),
            constant=True,
            sboTerm=SBO.KINETIC_CONSTANT,
            name="rate of a stage of the gastric emptying of the tablet",
            notes="""
            The gastric emptying of a tablet in the fasted state passes
            `N_EMPTYING` = 3 stages with this rate (fitted), the mean residence
            time in the stomach is 3 stages divided by the rate.
            """,
        ),
        Parameter(
            "Vmax_dis_riv",
            **FITS.parameter("GU__Vmax_dis_riv", unit=U.mg_per_hr),
            constant=True,
            sboTerm=SBO.MAXIMAL_VELOCITY,
            name="capacity of the dissolution rivaroxaban",
            notes="""
            The maximal rate of the dissolution of undissolved rivaroxaban in
            the intestine in the fasted state (solubility). The Cmax of the
            tablets rises less than the dose above 10 mg, 141, 173, 234 and
            316 µg/l at 10, 20, 40 and 80 mg (Kubitza2005a), the tmax of 10 mg
            is 2.75 hr (Kubitza2006).
            """,
            port=True,
        ),
        Parameter(
            "Km_dis_riv",
            1.0,
            U.mg,
            constant=True,
            sboTerm=SBO.MICHAELIS_CONSTANT,
            name="undissolved amount of half-maximal dissolution rivaroxaban",
            notes="""
            The dissolution saturates at a small undissolved amount: a dose of
            a tablet dissolves with the capacity `Vmax_dis_riv`, the remainder
            below `Km_dis_riv` with the first order rate `Vmax_dis_riv/Km_dis_riv`.
            """,
        ),
        Parameter(
            "ktransit_riv",
            **FITS.parameter("GU__ktransit_riv", unit=U.per_hr),
            constant=True,
            sboTerm=SBO.KINETIC_CONSTANT,
            name="rate of the transit of undissolved rivaroxaban",
            notes="""
            The undissolved rivaroxaban leaves the absorptive part of the
            intestine and is excreted with the feces, the small intestinal
            transit of about 3-4 hr; it limits the fraction of a high dose
            which dissolves.
            """,
        ),
        Parameter(
            "Mr_riv",
            435.88,
            U.g_per_mole,
            constant=True,
            name="Molecular weight rivaroxaban [g/mole]",
            sboTerm=SBO.MOLECULAR_MASS,
            port=True,
        ),
    ]
)

_m.rules.append(
    AssignmentRule(
        "application_riv",
        "Ka_application_riv/60 min_per_hr * PODOSE_riv/Mr_riv",
        unit=U.mmole_per_min,
    ),
)
_m.parameters.extend(
    [
        Parameter(
            f"emptying{k}_riv",
            np.nan,
            U.mmole_per_min,
            constant=False,
            sboTerm=SBO.RATE_OF_REACTION,
            name=f"gastric emptying rivaroxaban, stage {k}",
        )
        for k in range(1, N_EMPTYING + 1)
    ]
)
_m.rules.extend(
    [
        *(
            AssignmentRule(
                f"emptying{k}_riv",
                f"Ka_emptying_riv/60 min_per_hr * f_emptying * POSTOMACH{k}_riv/Mr_riv",
                unit=U.mmole_per_min,
            )
            for k in range(1, N_EMPTYING + 1)
        ),
        AssignmentRule(
            "POSTOMACH_riv",
            " + ".join(f"POSTOMACH{k}_riv" for k in range(1, N_EMPTYING + 1)),
            unit=U.mg,
        ),
    ]
)
_m.reactions.extend(
    [
        Reaction(
            sid="dissolution_riv",
            name="dissolution rivaroxaban",
            formula=(
                "f_solubility * Vmax_dis_riv/60 min_per_hr/Mr_riv"
                " * POINTESTINE_riv/(POINTESTINE_riv + Km_dis_riv)",
                U.mmole_per_min,
            ),
            equation="-> riv_lumen",
            compartment="Vlumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="""
            Dissolution of the undissolved dose in the intestine into the
            lumen; the rate saturates with the undissolved amount.
            """,
        ),
        Reaction(
            sid="transit_riv",
            name="transit undissolved rivaroxaban (feces)",
            formula=(
                "ktransit_riv/60 min_per_hr * POINTESTINE_riv/Mr_riv",
                U.mmole_per_min,
            ),
            equation="-> riv_feces",
            compartment="Vlumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="""
            The undissolved rivaroxaban leaves the absorptive part of the
            intestine and is excreted with the feces.
            """,
        ),
        Reaction(
            sid="application_solution_riv",
            name="application of the solution rivaroxaban",
            formula=(
                "Ka_application_riv/60 min_per_hr * POSOLUTION_riv/Mr_riv",
                U.mmole_per_min,
            ),
            equation="-> riv_lumen",
            compartment="Vlumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            notes="""
            The oral dose given as a solution enters the lumen of the intestine
            without dissolution.
            """,
        ),
    ]
)
_m.rate_rules.extend(
    [
        RateRule("PODOSE_riv", "-application_riv * Mr_riv", U.mg_per_min),
        RateRule(
            "POSOLUTION_riv",
            "-application_solution_riv * Mr_riv",
            U.mg_per_min,
        ),
        RateRule(
            "POSTOMACH1_riv",
            "(application_riv - emptying1_riv) * Mr_riv",
            U.mg_per_min,
        ),
        *(
            RateRule(
                f"POSTOMACH{k}_riv",
                f"(emptying{k - 1}_riv - emptying{k}_riv) * Mr_riv",
                U.mg_per_min,
            )
            for k in range(2, N_EMPTYING + 1)
        ),
        RateRule(
            "POINTESTINE_riv",
            f"(emptying{N_EMPTYING}_riv - dissolution_riv - transit_riv) * Mr_riv",
            U.mg_per_min,
        ),
    ]
)

model_intestine = _m


if __name__ == "__main__":
    from sbmlode import OdeSystem

    from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_system = OdeSystem.from_sbml(results.sbml_path)
    ode_system.write(results.sbml_path.parent / f"{results.sbml_path.stem}.md")
