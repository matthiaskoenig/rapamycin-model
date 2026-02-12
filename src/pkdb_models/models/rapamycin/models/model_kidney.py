"""Kidney model rapamycin."""

import numpy as np
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.rapamycin.models import annotations
from pkdb_models.models.rapamycin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "rapamycin_kidney"
version = 1

_m = Model(
    sid=mid,
    name="Model for renal rapamycin excretion.",
    notes=f"""
    Model for renal rapamycin  excretion.
    
    - ~2.2% of dose mainly as metabolites [Leung2006]

    **version** {version}
    
    ## Changelog
    
    **version 1**
    
    - initial model
        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),  # kidney
        (BQB.OCCURS_IN, "NCIT:C12415"),  # kidney

        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
        # (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
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
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
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
        "rx_urine",
        name="rapamycin metabolites (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True
    ),
]

_m.parameters.append(
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: reduced renal function
        """
    )
)

_m.reactions = [
    Reaction(
        sid="RXEX",
        name="rapamycin metabolites excretion (RXEX)",
        equation="rx_ext -> rx_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEX_k",
                1.0,
                U.per_min,
                name="rate urinary excretion of rapamycin metabolites",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * RXEX_k * Vki * rx_ext"
        )
    ),
]

model_kidney = _m


# def layout(dx=200, dy=200) -> pd.DataFrame:
#     """Layout definition."""
#
#     delta_x = 0.5 * dx
#     delta_y = 0.4 * dy
#
#     positions = [
#         ["ali_ext", 0, 0],
#         ["ALIEX", 0, 1.2 * delta_y],
#         ["ali_urine", 0, 2 * delta_y],
#     ]
#
#     df = pd.DataFrame(positions, columns=["id", "x", "y"])
#     df.set_index("id", inplace=True)
#     return df
#
# def annotations(dx=200, dy=200) -> list:
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
#         # kidney
#         cyviz.AnnotationShape(
#             x_pos=-delta_x, y_pos=-0.5 * delta_y, width=2 * delta_x, height=1* delta_y,
#             fill_color="#FF0000", **kwargs
#         ),
#         #plasma
#         cyviz.AnnotationShape(
#             x_pos=-delta_x, y_pos=0.5 * delta_y, width=2 * delta_x, height=1 * delta_y,
#             fill_color="#FFFFFF", **kwargs
#         ),
#         #urine
#         cyviz.AnnotationShape(
#             x_pos=-delta_x, y_pos=1.5 * delta_y, width=2 * delta_x, height=1 * delta_y,
#             fill_color="#000000", **kwargs
#         ),
#     ]
#     return annotations


if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    # SBML model
    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    # visualization
    from sbmlutils import cytoscape as cyviz
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    # cyviz.apply_layout(layout=layout())
    # cyviz.add_annotations(annotations=annotations())
    # cyviz.export_image(
    #     MODEL_BASE_PATH / f"{model_kidney.sid}.png",
    #     fit_content=True,
    # )
