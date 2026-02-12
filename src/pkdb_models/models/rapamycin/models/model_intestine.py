"""Rapamycin intestine model."""
import numpy as np
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.rapamycin.models import annotations
from pkdb_models.models.rapamycin.models import templates


class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "rapamycin_intestine",
    name="Model for rapamycin absorption in the small intestine",
    notes="""
    # Model for rapamycin absorption

    - large amount of radioactivity in feces: ~91.0% (mainly rapamycin metabolites)
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=annotations.model + [
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:45615"),  # gut
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "NCIT:C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine

        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        # (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ]
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
        "Vgu",
        1.2825,  # 0.0171 [l/kg] * 75 kg
        name="intestine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["gu"],
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
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Ventero",
        1.2825 * 0.1,  # 0.0171 [l/kg] * 75 kg * 0.9,
        name="intestinal lining (enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["apical"],
        spatialDimensions=2,
    ),
    Compartment(
        "Vbaso",
        np.nan,
        name="basolateral membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        annotations=annotations.compartments["basolateral"],
        spatialDimensions=2,
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
]

n_chain = 5

_m.parameters = [
    Parameter(
        "Vchain",
        0.1,
        unit=U.liter,
        name="volume of chain compartment",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),

]

for k in range(n_chain):
    # compartments for chain
    _m.compartments.append(
        Compartment(
            f"Vint_{k}",
            "Vchain",
            name="intestinal lumen (inner part of intestine)",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            annotations=annotations.compartments["gu_lumen"],
            constant=False
        ),
    )


_m.species = [
    Species(
        f"rap_stomach",
        metaId=f"meta_ali_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"rapamycin (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        boundaryCondition=True,
    ),
    Species(
        "rap_lumen",
        initialConcentration=0.0,
        name="rapamycin (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True,
    ),
    Species(
        "rx_lumen",
        initialConcentration=0.0,
        name="rapamycin metabolites (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
    Species(
        "rap_entero",
        initialConcentration=0.0,
        name="rapamycin (enterocytes)",
        compartment="Ventero",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True,
    ),
    Species(
        "rx_entero",
        initialConcentration=0.0,
        name="rapamycin metabolite (enterocytes)",
        compartment="Ventero",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
    Species(
        "rap_ext",
        initialConcentration=0.0,
        name="rapamycin (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True,
    ),
    Species(
        "rx_feces",
        initialConcentration=0.0,
        name="rapamycin metabolites (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),

]

for k in range(n_chain):
    # species for chain
    _m.species.append(
        Species(
            f"rx_int_{k}",
            initialConcentration=0.0,
            name=f"rapamycin metabolites (intestine) {k}",
            compartment=f"Vint_{k}",
            substanceUnit=U.mmole,
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["rx"],
        ),
    )


_m.reactions = [
    Reaction(
        "RAPIM",
        name="RAPIM",
        equation="rap_lumen -> rap_entero",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                "f_oatp",
                1.0,
                unit=U.dimensionless,
                name="OATP transporter activity",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""Used for drug-drug interactions.
                f_oatp = 1: normal activity;
                f_oatp < 1: reduced activity; inhibitor of absorption
                f_oatp > 1: increased activity; activator of absorption
                """
            ),
            Parameter(
                "RAPIM_k",
                0.10,
                unit=U.per_min,
                name="rate of rapamycin import enterocytes",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("f_oatp * RAPIM_k * Vgu * rap_lumen", U.mmole_per_min),
    ),
    Reaction(
        sid="RAPABS",
        name=f"absorption rapamycin (plasma)",
        compartment="Ventero",
        equation=f"rap_entero -> rap_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"RAPIM_k * Ventero * rap_entero",
            U.mmole_per_min),
    ),
    Reaction(
        sid="RAP2RX",
        name=f"rapamycin metabolism via CYP3A4/5",
        compartment="Ventero",
        equation=f"rap_entero -> rx_entero",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RAP2RX_k",
                0.02,
                unit=U.per_min,
                name="rate of rapamycin metabolism",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
            Parameter(
                "f_cyp3a4",
                1.0,
                unit=U.dimensionless,
                name="CYP3A4 activity",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""Used for drug-drug interactions.
                f_cyp3a4 = 1: normal activity;
                f_cyp3a4 < 1: reduced activity; inhibitor of CYP3A4
                f_cyp3a4 > 1: increased activity; activator of CYP3A4
                """
            ),
            Parameter(
                "f_cyp3a5",
                1.0,
                unit=U.dimensionless,
                name="CYP3A5 activity",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""Used for drug-drug interactions.
                f_cyp3a5 = 1: normal activity;
                f_cyp3a5 < 1: reduced activity; inhibitor of CYP3A4
                f_cyp3a5 > 1: increased activity; activator of CYP3A4
                """
            ),
        ],
        formula=(
            f"f_cyp3a4 * f_cyp3a5 * RAP2RX_k * Ventero * rap_entero",
            U.mmole_per_min),
    ),
    Reaction(
        "RXPG",
        name="PG",
        equation="rx_entero -> rx_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        pars=[
            Parameter(
                "f_pg",
                1.0,
                unit=U.dimensionless,
                name="P-gp transporter activity",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""Used for drug-drug interactions.
               f_pg = 1: normal activity;
               f_pg < 1: reduced activity; activator of absorption
               f_pg > 1: increased activity; inhibitor of absorption
               """
            ),
            Parameter(
                "RXPG_k",
                0.10,
                unit=U.per_min,
                name="rate of rapamycin metabolites PG export",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("f_pg * RXPG_k * Ventero * rx_entero", U.mmole_per_min),
    ),
    Reaction(
        sid="RXEXC",
        name=f"excretion rapamycin metabolites (feces)",
        compartment="Vlumen",
        equation=f"rx_lumen -> rx_int_0",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEXC_k",
                0.10,
                unit=U.per_min,
                name="rate of rapamycin metabolite fecal excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            f"RXEXC_k * Vgu  * rx_lumen",
            U.mmole_per_min),
    ),
]

for k in range(n_chain):
    # reactions for chain
    source = f"rx_int_{k}"
    target = f"rx_int_{k+1}" if k < n_chain - 1 else f"rx_feces"
    _m.reactions.append(
        Reaction(
            sid=f"RXEXC_{k}",
            name=f"excretion rapamycin metabolites {k}",
            compartment="Vlumen",
            equation=f"{source} -> {target}",
            sboTerm=SBO.TRANSPORT_REACTION,
            formula=(
                f"RXEXC_k * Vint_{k} * {source}",
                U.mmole_per_min,
            ),
        ),
    )


_m.parameters.extend([
    Parameter(
        f"PODOSE_rap",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose rapamycin [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_rap",
        2.0,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution rapamycin",
        port=True
    ),
    Parameter(
        f"Mr_rap",
        914.1719,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight rapamycin [g/mole]",
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
            sid=f"dissolution_rap",
            name=f"dissolution rapamycin",
            formula=(
                f"Ka_dis_rap/60 min_per_hr * PODOSE_rap/Mr_rap",
                U.mmole_per_min,
            ),
            equation=f"rap_stomach -> rap_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_rap", f"-dissolution_rap * Mr_rap", U.mg_per_min),
)


model_intestine = _m

# def layout(dx=200, dy=200) -> pd.DataFrame:
#     """Layout definition."""
#
#     delta_x = 0.5 * dx
#     delta_y = 0.4 * dy
#
#     positions = [
#         ["ali_stomach", 0, 0],
#         ["dissolution_ali", 0 * delta_x, 1 * delta_y],
#         ["ali_lumen", 0 * delta_x, 3 * delta_y],
#         ["ALIPG", 2 * delta_x, 2 * delta_y],
#         ["ALIIM", 2 * delta_x, 4 * delta_y],
#         ["ali_entero", 4 * delta_x, 3 * delta_y],
#         ["ALIMET", 5 * delta_x, 2 * delta_y],
#         ["ALIABS", 5 * delta_x, 4 * delta_y],
#         ["alm_entero", 7 * delta_x, 2 * delta_y],
#         ["ali_ext", 7 * delta_x, 4 * delta_y],
#
#         [f"ALIEXC", 0 * delta_x, 4 * delta_y],
#         ["ali_int_0", 0 * delta_x, 5 * delta_y],
#         ["ALIEXC_0", 0 * delta_x, 6 * delta_y],
#         ["ali_int_1", 1 * delta_x, 6 * delta_y],
#         ["ALIEXC_1", 2 * delta_x, 6 * delta_y],
#         ["ali_int_2", 3 * delta_x, 6 * delta_y],
#         ["ALIEXC_2", 3 * delta_x, 7 * delta_y],
#         ["ali_int_3", 3 * delta_x, 8 * delta_y],
#         ["ALIEXC_3", 2 * delta_x, 8 * delta_y],
#         ["ali_int_4", 1 * delta_x, 8 * delta_y],
#         ["ALIEXC_4", 0 * delta_x, 8 * delta_y],
#         ["ali_feces", 0 * delta_x, 9 * delta_y],
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
#         # intestine lumen (PG and OATP are at the apical side (it is near to the enterocytes)
#         cyviz.AnnotationShape(
#             x_pos=-1* delta_x, y_pos=- 0.5 * delta_y, width= 3 * delta_x, height=5 * delta_y,
#             fill_color="#FFFFFF", **kwargs
#         ),
#         #enterocytes
#         cyviz.AnnotationShape(
#             x_pos=3.5 * delta_x, y_pos= 1.25 * delta_y, width=4.5 * delta_x, height=2.25 * delta_y,
#             fill_color="#00FF00", **kwargs
#         ),
#         #plasma
#         cyviz.AnnotationShape(
#             x_pos=6 * delta_x, y_pos=3.5 * delta_y, width=2* delta_x, height=1 * delta_y,
#             fill_color="#FF0000", **kwargs
#         ),
#         # chain
#         cyviz.AnnotationShape(
#             x_pos=-1 * delta_x, y_pos=4.5 * delta_y, width=5* delta_x, height=4 * delta_y,
#             fill_color="#0000FF", **kwargs
#         ),
#         #feces
#         cyviz.AnnotationShape(
#             x_pos=-1 * delta_x, y_pos=8.5 * delta_y, width=2 * delta_x, height=1 * delta_y,
#             fill_color="#000000", **kwargs
#         ),
#
#     ]
#     return annotations


if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    # SBML model
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    # visualization
    from sbmlutils import cytoscape as cyviz
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
    # cyviz.apply_layout(layout=layout())
    # cyviz.add_annotations(annotations=annotations())
    # cyviz.export_image(
    #     MODEL_BASE_PATH / f"{model_intestine.sid}.png",
    #     fit_content=True,
    # )
