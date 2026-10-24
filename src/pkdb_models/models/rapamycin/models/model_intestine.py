"""Rapamycin intestine model."""

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

from pkdb_models.models.rapamycin.models import annotations, templates
from pkdb_models.models.rapamycin.models.fits import FITS


class U(templates.U):
    """Units of the intestine model."""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "rapamycin_intestine",
    name="Model for rapamycin absorption in the small intestine",
    notes="""
    # Model for rapamycin absorption

    The dose dissolves in the stomach and passes into the lumen of the
    intestine (`dissolution_rap`). A fraction `F_rap_abs` of the rapamycin
    which leaves the lumen is absorbed into the enterocytes, the rest is
    excreted with the feces (`RAPEXC`). In the enterocytes rapamycin is
    metabolized by CYP3A4 and CYP3A5 (`RAP2RX`), pumped back into the lumen by
    P-glycoprotein (`RAPPG`) or transported into the blood (`RAPBASO`); the
    metabolites formed in the enterocytes are released into the lumen and
    excreted with the feces after the transit through the colon.

    The fraction absorbed and the metabolism of the gut wall set the oral
    bioavailability of about 14 % (oral solution, Rapamune label); the
    radioactivity of an oral dose of 14C-rapamycin is recovered to 91 % in the
    feces, mainly as metabolites (Leung2006).

    Food (`fed`) slows the dissolution and the gastric emptying
    (`FED_f_dissolution`) and increases the fraction absorbed
    (`FED_f_absorption`): a high fat meal increases the AUC of rapamycin by
    35 % and delays its tmax from 0.8 to 3.1 hr (Zimmerman1999).
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
        name="blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["blood"],
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
        1.2825 * 0.1,  # 0.0171 [l/kg] * 75 kg * 0.1,
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

#: the compartments of the transit of the metabolites through the colon
n_chain = 5

_m.parameters = [
    Parameter(
        "f_entero",
        0.1,
        unit=U.dimensionless,
        name="fraction of the enterocytes of the intestine volume",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        "Vchain",
        0.1,
        unit=U.liter,
        name="volume of chain compartment",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
]
_m.rules = [
    AssignmentRule(
        "Ventero",
        "f_entero * Vgu",
        unit=U.liter,
        name="volume of the enterocytes",
    ),
]

for k in range(n_chain):
    # compartments for chain
    _m.compartments.append(
        Compartment(
            f"Vint_{k}",
            "Vchain",
            name=f"intestinal lumen {k} (inner part of intestine)",
            sboTerm=SBO.PHYSICAL_COMPARTMENT,
            unit=U.liter,
            annotations=annotations.compartments["gu_lumen"],
            constant=False,
        ),
    )


_m.species = [
    Species(
        "rap_stomach",
        metaId="meta_rap_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name="rapamycin (stomach)",
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
    ),
    Species(
        "rx_entero",
        initialConcentration=0.0,
        name="rapamycin metabolites (enterocytes)",
        compartment="Ventero",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rx"],
    ),
    Species(
        "rap_ext",
        initialConcentration=0.0,
        name="rapamycin (blood)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
        port=True,
    ),
    Species(
        "rap_colon",
        initialConcentration=0.0,
        name="rapamycin (colon)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["rap"],
    ),
    Species(
        "rap_feces",
        initialConcentration=0.0,
        name="rapamycin (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
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

# -------------------------------------
# Absorption
# -------------------------------------
_m.parameters.extend(
    [
        Parameter(
            "F_rap_abs",
            0.5,
            unit=U.dimensionless,
            name="fraction absorbed rapamycin (fasted)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Fraction of the rapamycin leaving the lumen which is absorbed into
            the enterocytes, the rest is excreted with the feces. Together with
            the metabolism of the gut wall it sets the oral bioavailability of
            about 14 % (Rapamune label); the split between the fraction
            absorbed and the gut wall cannot be identified from oral data and
            is assumed.
            """,
        ),
        Parameter(
            "RAPABS_k",
            0.05,
            unit=U.per_min,
            name="rate of rapamycin absorption",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "fed",
            0.0,
            unit=U.dimensionless,
            name="fed state (0 fasted, 1 fed)",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            The dose is taken in the fasted state (0) or with a high fat meal
            (1). The effects of food on the dissolution and gastric emptying
            and on the fraction absorbed are `FED_f_dissolution` and
            `FED_f_absorption`.
            """,
        ),
        Parameter(
            "FED_f_absorption",
            **FITS.parameter("GU__FED_f_absorption", unit=U.dimensionless),
            name="fed effect on the fraction absorbed",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FED_f_dissolution",
            **FITS.parameter("GU__FED_f_dissolution", unit=U.dimensionless),
            name="fed effect on the dissolution and gastric emptying",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_absorption",
            np.nan,
            unit=U.dimensionless,
            constant=False,
            name="factor of the fraction absorbed",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_dissolution",
            np.nan,
            unit=U.dimensionless,
            constant=False,
            name="factor of the dissolution and gastric emptying",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "absorption",
            np.nan,
            unit=U.mmole_per_min,
            constant=False,
            name="rapamycin leaving the lumen",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    ]
)
_m.rules.extend(
    [
        AssignmentRule(
            "f_absorption",
            "1 dimensionless + fed * (FED_f_absorption - 1 dimensionless)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "f_dissolution",
            "1 dimensionless + fed * (FED_f_dissolution - 1 dimensionless)",
            unit=U.dimensionless,
        ),
        AssignmentRule(
            "absorption",
            "RAPABS_k * Vlumen * rap_lumen",
            unit=U.mmole_per_min,
        ),
    ]
)

# -------------------------------------
# Gut wall
# -------------------------------------
_m.parameters.extend(
    [
        Parameter(
            "RAPBASO_k",
            0.05,
            unit=U.per_min,
            name="rate of rapamycin transport from the enterocytes into the blood",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "RAP2RX_k",
            0.1,
            unit=U.per_min,
            name="rate of rapamycin metabolism by CYP3A4 (gut wall)",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            Intrinsic clearance of CYP3A4 per volume of the enterocytes. The
            metabolism is linear: the Km of the microsomes (8.07 µM,
            Sattler1992) is above the concentrations of the enterocytes after
            the therapeutic doses.
            """,
        ),
        Parameter(
            "r_cyp3a5",
            1.0,
            unit=U.dimensionless,
            name="capacity of CYP3A5 relative to CYP3A4",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            port=True,
            notes="""
            Metabolism of rapamycin by CYP3A5 of the genotype *1/*1 relative
            to the one by CYP3A4.
            """,
        ),
        Parameter(
            "f_cyp3a4",
            1.0,
            unit=U.dimensionless,
            name="CYP3A4 activity",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Activity of CYP3A4 relative to the reference (*1/*1), changed by
            the genotype and by inhibitors and inducers:
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
            functional CYP3A5, the reference, most Caucasians), *1/*3 = 0.5,
            *1/*1 = 1.
            """,
        ),
        Parameter(
            "RAPPG_k",
            0.02,
            unit=U.per_min,
            name="rate of rapamycin efflux by P-glycoprotein",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "f_pg",
            1.0,
            unit=U.dimensionless,
            name="P-glycoprotein activity",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Activity of P-glycoprotein relative to the reference:
            f_pg = 1: normal activity;
            f_pg < 1: reduced activity, inhibition of the efflux;
            f_pg > 1: increased activity, induction of the efflux.
            """,
        ),
        Parameter(
            "RXEFF_k",
            1.0,
            unit=U.per_min,
            name="rate of the release of the metabolites into the lumen",
            sboTerm=SBO.KINETIC_CONSTANT,
        ),
        Parameter(
            "RXEXC_k",
            0.0033,
            unit=U.per_min,
            name="rate of the transit of the metabolites to the feces",
            sboTerm=SBO.KINETIC_CONSTANT,
            notes="""
            The metabolites pass the lumen and the chain of the colon in
            (n_chain + 1)/RXEXC_k, about 30 hr; the fecal radioactivity of
            14C-rapamycin is excreted within 24-72 hr (Leung2006).
            """,
        ),
    ]
)

_m.reactions = [
    Reaction(
        "RAPABS",
        name="absorption rapamycin (enterocytes)",
        equation="rap_lumen -> rap_entero",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_rap_abs * f_absorption * absorption", U.mmole_per_min),
    ),
    Reaction(
        "RAPEXC",
        name="excretion rapamycin (colon)",
        equation="rap_lumen -> rap_colon",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        formula=(
            "(1 dimensionless - F_rap_abs * f_absorption) * absorption",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "RAPFEC",
        name="excretion rapamycin (feces)",
        equation="rap_colon -> rap_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vlumen",
        # the mean transit of the metabolites, (n_chain + 1)/RXEXC_k
        formula=(
            f"RXEXC_k/{n_chain + 1} dimensionless * Vlumen * rap_colon",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        sid="RAPBASO",
        name="transport rapamycin (blood)",
        compartment="Vbaso",
        equation="rap_entero -> rap_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=("RAPBASO_k * Ventero * rap_entero", U.mmole_per_min),
    ),
    Reaction(
        sid="RAP2RX",
        name="rapamycin metabolism via CYP3A4/5",
        compartment="Ventero",
        equation="rap_entero -> rx_entero",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        formula=(
            "(f_cyp3a4 + r_cyp3a5 * f_cyp3a5) * RAP2RX_k * Ventero * rap_entero",
            U.mmole_per_min,
        ),
    ),
    Reaction(
        "RAPPG",
        name="efflux rapamycin by P-glycoprotein",
        equation="rap_entero -> rap_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("f_pg * RAPPG_k * Ventero * rap_entero", U.mmole_per_min),
    ),
    Reaction(
        "RXEFF",
        name="release rapamycin metabolites (lumen)",
        equation="rx_entero -> rx_lumen",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("RXEFF_k * Ventero * rx_entero", U.mmole_per_min),
    ),
    Reaction(
        sid="RXEXC",
        name="excretion rapamycin metabolites (feces)",
        compartment="Vlumen",
        equation="rx_lumen -> rx_int_0",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=("RXEXC_k * Vlumen * rx_lumen", U.mmole_per_min),
    ),
]

for k in range(n_chain):
    # reactions for chain
    source = f"rx_int_{k}"
    target = f"rx_int_{k + 1}" if k < n_chain - 1 else "rx_feces"
    _m.reactions.append(
        Reaction(
            sid=f"RXEXC_{k}",
            name=f"excretion rapamycin metabolites {k}",
            compartment=f"Vint_{k}",
            equation=f"{source} -> {target}",
            sboTerm=SBO.TRANSPORT_REACTION,
            formula=(
                f"RXEXC_k * Vint_{k} * {source}",
                U.mmole_per_min,
            ),
        ),
    )


_m.parameters.extend(
    [
        Parameter(
            "PODOSE_rap",
            0,
            U.mg,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="oral dose rapamycin [mg]",
            port=True,
        ),
        Parameter(
            "Ka_dis_rap",
            **FITS.parameter("GU__Ka_dis_rap", unit=U.per_hr),
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            name="Ka_dis [1/hr] dissolution rapamycin",
            port=True,
        ),
        Parameter(
            "Mr_rap",
            914.1719,
            U.g_per_mole,
            constant=True,
            name="Molecular weight rapamycin [g/mole]",
            sboTerm=SBO.MOLECULAR_MASS,
            port=True,
        ),
    ]
)

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid="dissolution_rap",
            name="dissolution rapamycin",
            formula=(
                "f_dissolution * Ka_dis_rap/60 min_per_hr * PODOSE_rap/Mr_rap",
                U.mmole_per_min,
            ),
            equation="rap_stomach -> rap_lumen",
            compartment="Vgu",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """,
        ),
    ]
)
_m.rate_rules.append(
    RateRule("PODOSE_rap", "-dissolution_rap * Mr_rap", U.mg_per_min),
)


model_intestine = _m

if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    # SBML model
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )
