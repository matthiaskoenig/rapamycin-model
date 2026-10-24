"""Liver model rapamycin."""

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

from pkdb_models.models.rapamycin.models import annotations, templates
from pkdb_models.models.rapamycin.models.fits import FITS


class U(templates.U):
    """Units of the liver model."""


mid = "rapamycin_liver"

_m = Model(
    sid=mid,
    name="Model for hepatic rapamycin metabolism.",
    notes="""
    # Model for hepatic rapamycin metabolism

    Rapamycin equilibrates fast between the blood and the liver and is
    metabolized by CYP3A4 and CYP3A5 (`RAP2RX`), linear in the concentrations
    of the therapeutic doses. The concentrations are those of whole blood:
    about 95 % of the rapamycin of the blood is bound in the erythrocytes
    (FKBP12), the intrinsic clearance is a clearance of whole blood. Rapamycin
    has a low hepatic extraction (blood clearance about 3 l/hr against a
    hepatic blood flow of about 90 l/hr).

    The metabolites are released into the blood (`RXEX`) and excreted into the
    bile (`RXBEX`), from where they reach the lumen of the intestine and the
    feces.
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "ncit/C12392"),
        (BQB.HAS_PROPERTY, "ncit/C79371"),  # Pharmacokinetics: Metabolism
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
        name="rapamycin (blood)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True,
    ),
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
        0.15,
        unit=U.per_min,
        name="rate of the biliary excretion of the rapamycin metabolites",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""
        The biliary clearance of the metabolites is about 5 times their renal
        clearance (`KI__RXEX_k`): after an oral dose of 14C-rapamycin 91 % of
        the radioactivity is recovered in the feces and 2.2 % in the urine
        (Leung2006), and most of the fecal radioactivity is formed in the gut
        wall or not absorbed.
        """,
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
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("RAPIM_k * Vli * (rap_ext - rap)", U.mmole_per_min),
        notes="""Assumed to be fast.""",
    ),
    Reaction(
        sid="RAP2RX",
        name="rapamycin metabolism (CYP3A4/5)",
        equation="rap -> rx",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "RAP2RX_k",
                **FITS.parameter("LI__RAP2RX_k", unit=U.per_min),
                name="rate of rapamycin metabolism by CYP3A4",
                sboTerm=SBO.KINETIC_CONSTANT,
                notes="""
                Intrinsic clearance of CYP3A4 per volume of the liver, of
                whole blood concentrations. The metabolism is linear: the Km of
                the microsomes (8.07 µM, Sattler1992) is far above the
                concentrations of the therapeutic doses.
                """,
            ),
            Parameter(
                "r_cyp3a5",
                1.0,
                U.dimensionless,
                name="capacity of CYP3A5 relative to CYP3A4",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                port=True,
            ),
            Parameter(
                "f_cyp3a4",
                1.0,
                unit=U.dimensionless,
                name="CYP3A4 activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                Activity of CYP3A4 relative to the reference (*1/*1), changed
                by the genotype and by inhibitors and inducers:
                f_cyp3a4 = 1: normal activity;
                f_cyp3a4 < 1: reduced activity, inhibition of CYP3A4;
                f_cyp3a4 > 1: increased activity, induction of CYP3A4.
                """,
            ),
            Parameter(
                "f_cyp3a5",
                0.0,
                unit=U.dimensionless,
                name="CYP3A5 activity",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                The fraction of functional CYP3A5 alleles: *3/*3 = 0 (no
                functional CYP3A5, the reference, most Caucasians), *1/*3 =
                0.5, *1/*1 = 1.
                """,
            ),
        ],
        formula=(
            "(f_cyp3a4 + r_cyp3a5 * f_cyp3a5) * RAP2RX_k * Vli * rap",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="RXEX",
        name="rapamycin metabolites export (RXEX)",
        equation="rx <-> rx_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "RXEX_k",
                100.0,
                U.per_min,
                name="rate rapamycin metabolites export",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=("RXEX_k * Vli * (rx - rx_ext)", U.mmole_per_min),
        notes="""Assumed to be fast.""",
    ),
    Reaction(
        "RXBEX",
        name="rapamycin metabolites biliary excretion",
        equation="rx -> rx_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=(
            "RXBEX_k * Vli * rx",
            U.mmole_per_min,
        ),
        notes="""
        The metabolites are excreted into the bile, which reaches the lumen of
        the intestine; they are not reabsorbed.
        """,
    ),
]

model_liver = _m

if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )
