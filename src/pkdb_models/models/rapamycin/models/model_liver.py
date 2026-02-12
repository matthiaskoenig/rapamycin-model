"""Liver model rapamycin."""

import numpy as np
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.rapamycin.models import annotations
from pkdb_models.models.rapamycin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "rapamycin_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic rapamycin metabolism.",
    notes=f"""
    # Model for rapamycin metabolism.
    
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "NCIT:C12392"),

        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ]
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
        name="plasma membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma membrane"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
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
        1.2825 * 0.9,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lumen (inner part of intestine)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),

]

_m.species = [
    Species(
        "rap_ext",
        name="rapamycin (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True
    ),
    Species(
        "rx_ext",
        name="rapamycin metabolites (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True
    ),
    Species(
        "rap",
        name="rapamycin (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
    ),
    Species(
        "rx",
        name="rapamycin metabolites (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
    ),
    Species(
        "rx_bi",
        initialConcentration=0.0,
        name="rapamycin metabolites (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        notes="""
        Bile rapamycin metabolites in amount.
        """,
    ),
    Species(
        "rx_lumen",
        initialConcentration=0.0,
        name="rapamycin metabolites (lumen)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
]

_m.parameters = [
    Parameter(
        "RXBEX_k",
        0.0001,
        unit=U.per_min,
        name="rate for rapamycin metabolites export in bile",
        sboTerm=SBO.KINETIC_CONSTANT,
    )
]

_m.reactions = [
    Reaction(
        sid="RAPIM",
        name="rapamycin import (RAPIM)",
        equation="rap_ext <-> rap",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RAPIM_k",
                100.0,
                U.per_min,
                name="rate rapamycin import",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        ],
        formula=(
            "RAPIM_k * Vli * (rap_ext - rap)"
        ),
        notes = """Assumed to be fast."""
    ),
    Reaction(
        sid="RAP2RX",
        name="rapamycin metabolism (CYP3A4/5)",
        equation="rap -> rx",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "RAP2RX_Vmax",
                0.02,
                U.mmole_per_min_l,
                name="Vmax rapamycin metabolism",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
            Parameter(
                "RAP2RX_Km_rap",
                2.9E-3,     # 2.9 [μM]
                U.mM,
                name="Km rapamycin metabolism",
                sboTerm=SBO.MICHAELIS_CONSTANT,
                notes="""Michaelis-Menten constant for rapamycin metabolism from in vitro studies.
                
                Km value in the range of approximately 1.1μM to 4.7μM for sirolimus. [perplexity]
                """
            ),
            Parameter(
                "f_cyp3a4",
                1,
                U.dimensionless,
                name="scaling factor cyp3a4 activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""Scaling factor to vary CYP3A4 activity.
                1.0: unchanged activity; < 1.0 decreased activity; >1.0 increased activity.
                """
            ),
            Parameter(
                "f_cyp3a5",
                1,
                U.dimensionless,
                name="scaling factor cyp3a5 activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""Scaling factor to vary CYP3A5 activity.
                1.0: unchanged activity; < 1.0 decreased activity; >1.0 increased activity.
                """
            ),
        ],
        formula=(
            "f_cyp3a4 * f_cyp3a5 * RAP2RX_Vmax * Vli * rap/(rap + RAP2RX_Km_rap)"
        ),
    ),
    Reaction(
        sid="RXEX",
        name="rapamycin export (RXEX)",
        equation="rx <-> rx_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEX_k",
                100.0,
                U.per_min,
                name="rate rapamycin metabolites export",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            ),
        ],
        formula=(
            "RXEX_k * Vli * (rx - rx_ext)"
        ),
        notes = """Assumed to be fast."""
    ),

    Reaction(
        "RXBEX",
        name="rapamycin metabolites bile export",
        equation="rx -> rx_bi",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=(
            "RXBEX_k * Vli * rx",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "RXEHC",
        name="rapamycin metabolites enterohepatic circulation",
        equation="rx_bi -> rx_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        formula=("RXEX", U.mmole_per_min),
    ),
]

# def liver_layout(dx=200, dy=200) -> pd.DataFrame:
#     """Layout definition."""
#
#     delta_x = 0.5 * dx
#     delta_y = 0.4 * dy
#
#     positions = [
#         ["ali_ext", 0, 0],
#         ["ALIIM", 0, 1 * delta_y],
#         ["ali", 0, 2 * delta_y],
#         ["ALM", 0, 3 * delta_y],
#         ["alm", 0, 4 * delta_y],
#
#         ["ALIEX", delta_x, 2 * delta_y],
#         ["ali_bi", 2*delta_x, 2 * delta_y],
#         ["ALIEHC", 2*delta_x, 1 * delta_y],
#         ["ali_lumen", 2 * delta_x, 0 * delta_y],
#     ]
#
#     df = pd.DataFrame(positions, columns=["id", "x", "y"])
#     df.set_index("id", inplace=True)
#     return df
#
# def liver_annotations(dx=200, dy=200) -> list:
#
#     delta_x = 0.5 * dx
#     delta_y = 0.4 * dy
#
#     kwargs = {
#         "type": cyviz.AnnotationShapeType.ROUND_RECTANGLE,
#         "opacity": 20,
#         "border_color": "#000000",
#         "border_thickness": 2,
#     }
#     annotations = [
#         # liver
#         cyviz.AnnotationShape(
#             x_pos=-delta_x, y_pos=delta_y, width=2 * delta_x, height=3.5* delta_y,
#             fill_color="#FFFFFF", **kwargs
#         ),
#         #plasma
#         cyviz.AnnotationShape(
#             x_pos=-delta_x, y_pos=- 0.75 * delta_y, width=2 * delta_x, height=1.75 * delta_y,
#             fill_color="#FF0000", **kwargs
#         ),
#         #intestine
#         cyviz.AnnotationShape(
#             x_pos=delta_x, y_pos=- 0.75 * delta_y, width=2 * delta_x, height=1.75 * delta_y,
#             fill_color="#FFFFFF", **kwargs
#         ),
#         #bile duct
#         cyviz.AnnotationShape(
#             x_pos=delta_x, y_pos= 1 * delta_y, width=2 * delta_x, height=2 * delta_y,
#             fill_color="#000000", **kwargs
#         ),
#     ]
#     return annotations


model_liver = _m

if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    # visualization
    from sbmlutils import cytoscape as cyviz
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    # cyviz.apply_layout(layout=liver_layout())
    # cyviz.add_annotations(annotations=liver_annotations())
    # cyviz.export_image(
    #     MODEL_BASE_PATH / f"{model_liver.sid}.png",
    #     fit_content=True,
    # )