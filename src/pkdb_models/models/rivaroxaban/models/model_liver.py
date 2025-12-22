"""Liver model for rivaroxaban."""

from sbmlutils import cytoscape as cyviz
from pathlib import Path
import pandas as pd
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.rivaroxaban.models import annotations
from pkdb_models.models.rivaroxaban.models import templates
from pkdb_models.models.rivaroxaban import MODEL_BASE_PATH


class U(templates.U):
    """UnitDefinitions"""

    pass

# TODO: add EHC equations
mid = "rivaroxaban_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic rivaroxaban metabolism.",
    notes=f"""
    Model for hepatic rivaroxaban metabolism.
    
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "ncit/C12392"),

        (BQB.HAS_PROPERTY, "ncit/C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ] + templates.model_annotations
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
        port=True
    ),
    Compartment(
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True
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
        port=True
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
        port=True
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
        compartment="Vli",
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
        9.99999999999998e-05,
        unit=U.per_min,
        name="rate for rivaroxaban metabolite export in bile",
        sboTerm=SBO.KINETIC_CONSTANT,
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
                "RIVIM_Vmax",
                1000.0,
                U.per_min,
                name="Vmax rivaroxaban import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "RIVIM_Vmax * Vli * (riv_ext - riv)"
        ),
        notes="""Assuming fast equilibration."""
    ),
    Reaction(
        sid="RIV2RX",
        name="rivaroxaban conversion (RIV2RX)",
        equation="riv -> rx",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "RIV2RX_Vmax",
                0.06401213063939644,
                U.per_min,
                name="Vmax rivaroxaban conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "RIV2RX_Vmax * Vli * riv"
        ),
        notes="""
        Combination of multiple enzymes in the liver which create a variet of different
        metabolites. Not much is known about the conversions of rivaroxaban to its metabolites.
        [Weinz2009].
        """
    ),
    Reaction(
        sid="RXEX",
        name="rivaroxaban metabolite export (RXEX)",
        equation="rx <-> rx_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEX_Vmax",
                1000.0,
                U.per_min,
                name="Vmax rivaroxaban metabolite export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "RXEX_Vmax * Vli * (rx - rx_ext)"
        ),
        notes="""Assuming fast equilibration."""
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

def liver_layout(dx=200, dy=200) -> pd.DataFrame:
    """Liver layout for Cytoscape visualization."""
    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    positions = [
        ["riv_ext", delta_x, 0.5 * delta_y],
        ["RIVIM", delta_x, 1.5 * delta_y],
        ["riv", 1.2 * delta_x, 2.2 * delta_y],
        ["RIV2RX", 1.6 * delta_x, 3.1 * delta_y],
        ["rx", 2 * delta_x, 4 * delta_y],

        ["RXEX", 3 * delta_x, 1.5 * delta_y],
        ["rx_ext", 3 * delta_x, 0.5 * delta_y],

        ["RXEX_BI", 2 * delta_x, 5.5 * delta_y],
        ["rx_bi", 2 * delta_x, 6.3 * delta_y],

        ["RXEXEHC", 2 * delta_x, 7.5 * delta_y],
        ["rx_lumen", 2 * delta_x, 8.5 * delta_y],
    ]

    df = pd.DataFrame(positions, columns=["id", "x", "y"])
    df.set_index("id", inplace=True)
    return df


def liver_annotations(dx=200, dy=200) -> list:
    """Bounding boxes for 'plasma', 'liver', 'bile', and 'gut lumen' compartments stacked vertically."""

    delta_x = 0.8 * dx
    delta_y = 0.4 * dy

    width = 4 * delta_x
    x_pos = 0

    kwargs = {
        "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
        "opacity": 20,
        "border_color": "#000000",
        "border_thickness": 2,
    }

    annotations = [
        # Plasma
        cyviz.AnnotationShape(
            x_pos=x_pos,
            y_pos=0 * delta_y,
            width=width,
            height=1.5 * delta_y,
            fill_color="#FF0000",  # red
            **kwargs
        ),

        # Liver
        cyviz.AnnotationShape(
            x_pos=x_pos,
            y_pos=1.5 * delta_y,
            width=width,
            height=4 * delta_y,
            fill_color="#FFFFFF",  # white
            **kwargs
        ),

        # Bile
        cyviz.AnnotationShape(
            x_pos=x_pos,
            y_pos=5.5 * delta_y,
            width=width,
            height=2 * delta_y,
            fill_color="#00FF00",  # green
            **kwargs
        ),

        # Gut lumen
        cyviz.AnnotationShape(
            x_pos=x_pos,
            y_pos=7.5 * delta_y,
            width=width,
            height=1.5 * delta_y,
            fill_color="#FFE3E3",  # light pink
            **kwargs
        ),
    ]

    return annotations



if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    cyviz.apply_layout(layout=liver_layout())
    cyviz.add_annotations(annotations=liver_annotations())


# PNG output
#     figure_path: Path = MODEL_BASE_PATH.parent / "figures" / f"{model_liver.sid}.png"
#     cyviz.export_image(
#         figure_path,
#         fit_content=True,
#      )
#     console.print(figure_path)