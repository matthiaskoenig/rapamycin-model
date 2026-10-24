"""Kidney model rapamycin."""

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

from pkdb_models.models.rapamycin.models import annotations, templates


class U(templates.U):
    """Units of the kidney model."""


mid = "rapamycin_kidney"
_m = Model(
    sid=mid,
    name="Model for renal rapamycin excretion.",
    notes="""
    # Model for the renal excretion of the rapamycin metabolites

    Rapamycin is not excreted by the kidneys; its metabolites are excreted in
    the urine: 2.2 % of the radioactivity of an oral dose of 14C-rapamycin is
    recovered in the urine (Leung2006).
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),  # kidney
        (BQB.OCCURS_IN, "ncit/C12415"),  # kidney
        (BQB.HAS_PROPERTY, "ncit/C79372"),  # Pharmacokinetics: Excretion
    ],
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="blood",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["blood"],
        port=True,
    ),
    Compartment(
        "Vki",
        value=0.3,  # 0.4 % of bodyweight
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
        "rx_ext",
        name="rapamycin metabolites (blood)",
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
        name="rapamycin metabolites (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
        port=True,
    ),
]

_m.parameters.append(
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        notes="""scaling factor for renal function. 1.0: normal renal function;
        <1.0: reduced renal function
        """,
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
                0.12,
                U.per_min,
                name="rate urinary excretion of rapamycin metabolites",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("f_renal_function * RXEX_k * Vki * rx_ext", U.mmole_per_min),
    ),
]

model_kidney = _m


if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    # SBML model
    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )
