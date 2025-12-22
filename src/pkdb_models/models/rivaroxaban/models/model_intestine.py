"""Rivaroxaban intestine model."""

from pathlib import Path
import pandas as pd
from sbmlutils.factory import create_model, FactoryResult
from sbmlutils.converters import odefac
from sbmlutils import cytoscape as cyviz
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.rivaroxaban.models import annotations
from pkdb_models.models.rivaroxaban.models import templates
from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH

class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "rivaroxaban_intestine",
    name="Model for rivaroxaban absorption in the small intestine",
    notes="""
    # Model for rivaroxaban absorption
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:45615"),  # gut
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "ncit/C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine

        (BQB.HAS_PROPERTY, "ncit/C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ] + templates.model_annotations
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
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
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
        f"riv_stomach",
        metaId=f"meta_riv_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"rivaroxaban (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["riv"],
        boundaryCondition=True,
    ),
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

_m.parameters = [
    Parameter(
        f"F_riv_abs",
        0.82,  # [0.82 - 0.88]
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed rivaroxaban",
        notes="""
        Fraction absorbed, i.e., only a fraction of the rivaroxaban in the intestinal lumen
        is absorbed. This parameter determines how much of the rivaroxaban is excreted.

        `F_riv_abs` of dose is absorbed. `(1-F_riv_abs)` is excreted in feces.

        ~12 % in feces and 6% unaccounted, total around 18 %
        applicant’s estimate of 80% to 100% appears.
        
        assuming: 0.82 in fasted state
        assuming: 1.00 in fed state
        """,
    ),
    Parameter(
        "RIVABS_k",
        0.052027850297051675,
        unit=U.per_min,
        name="rate of rivaroxaban absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "f_absorption",
        1,
        unit=U.dimensionless,
        name="scaling factor for absorption rate",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""1.0: normal absorption corresponding to tablet under fasting conditions.

        allows to change the velocity of absorption
        """
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption_riv",
        value="f_absorption * RIVABS_k * Vlumen * riv_lumen",
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
        name=f"excretion rivaroxaban (feces)",
        compartment="Vlumen",
        equation=f"riv_lumen -> riv_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"(1 dimensionless - F_riv_abs) * absorption_riv",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="RXEXC",
        name=f"excretion rivaroxaban metabolites (feces)",
        compartment="Vlumen",
        equation=f"rx_lumen -> rx_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEXC_k",
                0.10,
                unit=U.per_min,
                name="rate of rivaroxaban metabolites excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            f"f_absorption * RXEXC_k * Vlumen * rx_lumen",
            U.mmole_per_min,
        ),
    ),
]

_m.parameters.extend([
    Parameter(
        f"PODOSE_riv",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose rivaroxaban [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_riv",
        0.15960999342060214,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution rivaroxaban",
        port=True
    ),
    Parameter(
        f"Mr_riv",
        435.88,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight rivaroxaban [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_riv",
            name=f"dissolution rivaroxaban",
            formula=(
                f"Ka_dis_riv/60 min_per_hr * PODOSE_riv/Mr_riv",
                U.mmole_per_min,
            ),
            equation=f"riv_stomach -> riv_lumen",
            compartment="Vlumen",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_riv", f"-dissolution_riv * Mr_riv", U.mg_per_min),
)

model_intestine = _m


def intestine_layout(dx=200, dy=200) -> pd.DataFrame:
    """Intestine layout for cytoscape visualization."""
    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    positions = [
        ["riv_stomach",      0 * delta_x, 0],
        ["dissolution_riv",  0 * delta_x, 1 * delta_y],
        ["riv_lumen",        0.6 * delta_x, 2 * delta_y],
        ["RIVABS",           0.6 * delta_x, 3 * delta_y],
        ["riv_ext",          0.6 * delta_x, 3.75 * delta_y],
        ["RIVEXC",           1.2 * delta_x, 1 * delta_y],
        ["riv_feces",        1.2 * delta_x, 0 * delta_y],
        ["rx_lumen",         2.3 * delta_x, 2 * delta_y],
        ["RXEXC",            2.3 * delta_x, 1 * delta_y],
        ["rx_feces",         2.3 * delta_x, 0 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df



def intestine_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'stomach', 'intestines', 'plasma' and 'feces'."""
    from sbmlutils.cytoscape import AnnotationShape, AnnotationShapeType

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    kwargs = dict(
        type=AnnotationShapeType.ROUND_RECTANGLE,
        opacity=20,
        border_color="#000000",
        border_thickness=2,
    )

    # Stomach region
    stomach_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=-0.5 * delta_y,
        width=1.2 * delta_x,
        height=1.5 * delta_y,
        fill_color="#FFC0CB",      # pinkish
        **kwargs
    )

    # Intestinal lumen
    lumen_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=stomach_box.y_pos + stomach_box.height,
        width=3.6 * delta_x,
        height=2 * delta_y,
        fill_color="#ffe3e3",      # very light pink
        **kwargs
    )

    # Plasma
    plasma_box = AnnotationShape(
        x_pos=-0.6 * delta_x,
        y_pos=lumen_box.y_pos + lumen_box.height,
        width=3.6 * delta_x,
        height=1.25 * delta_y,
        fill_color="#FF0000",  # pale pink
        **kwargs
    )

    # Feces
    feces_box = AnnotationShape(
        x_pos=0.6 * delta_x,
        y_pos=-0.5 * delta_y,
        width=2.4 * delta_x,
        height=1.5 * delta_y,
        fill_color="#d7c2af",  # beige
        **kwargs
    )

    return [stomach_box, lumen_box, plasma_box, feces_box]



if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=False)
    cyviz.apply_layout(layout=intestine_layout())
    cyviz.add_annotations(annotations=intestine_annotations())


    # PNG output
    # figure_path: Path = MODEL_BASE_PATH.parent / "figures" / f"{model_intestine.sid}.png"
    # cyviz.export_image(
    #     figure_path,
    #     fit_content=True,
    #  )
    # console.print(figure_path)