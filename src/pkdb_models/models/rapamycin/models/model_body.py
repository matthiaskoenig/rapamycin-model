"""Rapamycin body model."""

from typing import NotRequired, TypedDict

from sbmlutils.factory import (
    PORT_SUFFIX,
    AssignmentRule,
    Compartment,
    ExternalModelDefinition,
    Model,
    Parameter,
    RateRule,
    Reaction,
    ReplacedElement,
    Species,
    Submodel,
    UnitDefinition,
    create_model,
)
from sbmlutils.metadata import (
    BQB,
    SBO,
)

from pkdb_models.models.rapamycin.models import annotations, templates
from pkdb_models.models.rapamycin.models.fits import FITS


# -----------------------------------------------------------------------------
# Whole Body Metabolism
# -----------------------------------------------------------------------------
class U(templates.U):
    """UnitsDefinitions."""

    cm = UnitDefinition("cm")
    kg = UnitDefinition("kg")

    l_per_kg = UnitDefinition("l_per_kg", "l/kg")

    mg_per_min = UnitDefinition("mg_per_min", "mg/min")

    ml = UnitDefinition("ml")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_s = UnitDefinition("ml_per_s", "ml/s")
    ml_per_s_kg = UnitDefinition("ml_per_s_kg", "ml/s/kg")

    s_per_min = UnitDefinition("s_per_min", "s/min")


_m = Model(
    "rapamycin_body",
    name="Rapamycin body model",
    notes="""
    # Rapamycin body model

    Whole-body model of rapamycin (sirolimus) with the blood circulation, the
    tissue distribution and the submodels of the intestine, the liver and the
    kidneys.

    The circulating species are whole blood concentrations: rapamycin is
    bound in the erythrocytes (FKBP12, blood to plasma ratio 36-142), nearly
    all studies measure whole blood. The blood flows transport the blood
    concentrations, the plasma concentration is the blood concentration over
    the blood to plasma ratio `BP_rap`.
    """
    + templates.terms_of_use,
    units=U,
    model_units=templates.model_units,
    creators=templates.creators,
    annotations=[
        *annotations.model,
        # tissue
        (BQB.OCCURS_IN, "BTO:0001489"),  # whole body
        (BQB.OCCURS_IN, "FMA:256135"),  # body
        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
        (BQB.HAS_PROPERTY, "NCIT:C79370"),  # Pharmacokinetics: Distribution
    ],
)

_m.ports = []
_m.deletions = []
_m.replaced_elements = []


class SubstanceBody(TypedDict):
    """Definition of a substance transported in the whole-body model."""

    name: str
    unit: UnitDefinition
    cinit: float  # initial concentration [mmole/l]
    Mr: NotRequired[float]  # molecular weight [g/mole]
    IVDOSE: NotRequired[float]  # intravenous dose [mg]
    PODOSE: NotRequired[float]  # oral dose [mg]
    #: distribution into the tissues (`ftissue_<sid>`, `Kp_<sid>`) and slow
    #: binding in the tissue (`Kpdeep_<sid>`, `kdeep_<sid>`), fitted
    tissue_distribution: NotRequired[bool]
    BP: NotRequired[float]  # blood to plasma ratio [-]


SUBSTANCES_BODY: dict[str, SubstanceBody] = {
    "rap": {
        "name": "rapamycin",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # Molecular weight
        "Mr": 914.1719,  # [g/mole]
        # iv doses [mg]
        "IVDOSE": 0,  # venous blood
        # oral absorption
        "PODOSE": 0,  # dose
        # tissue distribution
        # the distribution in tissues (ftissue, Kp) and the slow binding in
        # the tissue (Kpdeep, kdeep) are fitted, see `fits.FITS`
        # logP: 4.63, https://go.drugbank.com/drugs/DB00877
        "tissue_distribution": True,
        # blood to plasma ratio: 142 ± 39 in healthy subjects (Leung2006),
        # 84.5-117 (Zimmerman2008), 36 ± 12 in transplant patients
        # (Brattstrom2000)
        "BP": 142,
    },
    "rx": {
        "name": "rapamycin metabolites",
        "unit": U.mmole,
        # initial concentration
        "cinit": 0.0,  # [mmole/l]
        # the metabolites are not dosed, i.e., no molecular weight, and are
        # not distributed into the tissues
    },
}

# -----------------------------------------------------------------------------
# Submodels
# -----------------------------------------------------------------------------
COMPARTMENTS_BODY = {
    "ki": "kidney",
    "li": "liver",
    "lu": "lung",
    "gu": "gut",
    "re": "rest",
    "ar": "arterial blood",
    "ve": "venous blood",
    "po": "portal vein",
    "hv": "hepatic vein",
}

SUBMODEL_SID_DICT = {  # tissue to submodel mapping
    "ki": "KI",  # kidney
    "li": "LI",  # liver
    "gu": "GU",  # gut/intestine
}

kidney_id = "rapamycin_kidney"
liver_id = "rapamycin_liver"
gut_id = "rapamycin_intestine"

_m.external_model_definitions = [
    ExternalModelDefinition(
        sid="kidney",
        name="kidney model definition",
        source=f"{kidney_id}.xml",
        modelRef=kidney_id,
    ),
    ExternalModelDefinition(
        sid="liver",
        name="liver model definition",
        source=f"{liver_id}.xml",
        modelRef=liver_id,
    ),
    ExternalModelDefinition(
        sid="gut", name="gut model definition", source=f"{gut_id}.xml", modelRef=gut_id
    ),
]

_m.submodels = [
    Submodel(sid=SUBMODEL_SID_DICT["ki"], modelRef="kidney", name="kidney submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["li"], modelRef="liver", name="liver submodel"),
    Submodel(sid=SUBMODEL_SID_DICT["gu"], modelRef="gut", name="gut submodel"),
]

# -------------------------------------------------------------------------------------------------
# Compartments
# -------------------------------------------------------------------------------------------------
_m.compartments = [
    Compartment(
        "Vre",
        metaId="meta_Vre",
        value=1,
        unit=U.liter,
        constant=False,
        name="rest of body",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["re"],
    ),
    Compartment(
        "Vgu",
        metaId="meta_Vgu",
        value=1,
        unit=U.liter,
        constant=False,
        name="gut",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu"],
    ),
    Compartment(
        "Vki",
        metaId="meta_Vki",
        value=1,
        unit=U.liter,
        constant=False,
        name="kidney",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
    ),
    Compartment(
        "Vli",
        metaId="meta_Vli",
        value=1,
        unit=U.liter,
        constant=False,
        name="liver",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
    ),
    Compartment(
        "Vlu",
        metaId="meta_Vlu",
        value=1,
        unit=U.liter,
        constant=False,
        name="lung",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["lu"],
    ),
    Compartment(
        "Vve",
        metaId="meta_Vve",
        value=1,
        unit=U.liter,
        constant=False,
        name="venous blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ve"],
    ),
    Compartment(
        "Var",
        metaId="meta_Var",
        value=1,
        unit=U.liter,
        constant=False,
        name="arterial blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ar"],
    ),
    Compartment(
        "Vurine",
        metaId="meta_Vurine",
        value=1,
        unit=U.liter,
        constant=True,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["urine"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
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
        annotations=annotations.compartments["stomach"],
    ),
    Compartment(
        "Vpo",
        metaId="meta_Vpo",
        value=1,
        unit=U.liter,
        constant=False,
        name="portal blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["po"],
    ),
    Compartment(
        "Vhv",
        metaId="meta_Vhv",
        value=1,
        unit=U.liter,
        constant=False,
        name="hepatic venous blood",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["hv"],
    ),
]

# create blood and tissue compartments (for correct blood volume)
for cid in COMPARTMENTS_BODY:
    if cid not in ["ar", "ve", "po", "hv"]:
        _m.compartments.extend(
            [
                Compartment(
                    f"V{cid}_tissue",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} tissue",
                    metaId=f"meta_V{cid}_tissue",
                    port=True,
                    annotations=annotations.compartments["parenchyma"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
                Compartment(
                    f"V{cid}_blood",
                    value=1,
                    unit=U.liter,
                    constant=False,
                    name=f"{cid} blood",
                    metaId=f"meta_V{cid}_blood",
                    port=True,
                    annotations=annotations.compartments["blood"]
                    + [
                        (BQB.IS_PART_OF, aid[1])
                        for aid in annotations.compartments[cid]
                    ],
                    sboTerm=SBO.PHYSICAL_COMPARTMENT,
                ),
            ]
        )

# calculate gut lumen
_m.compartments.append(
    Compartment(
        "Vgulumen",
        metaId="meta_Vgulumen",
        value=1,
        unit=U.liter,
        constant=False,
        name="gut lumen",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
    )
)
_m.parameters = [
    Parameter("f_lumen", 0.9, name="fraction lumen of intestine", unit=U.dimensionless)
]
_m.rules = [
    AssignmentRule("Vgulumen", "f_lumen * Vgu", unit=U.liter),
]


# replace volumes of submodels
_m.replaced_elements.extend(
    [
        # kidney
        ReplacedElement(
            sid="Vki_tissue_RE",
            metaId="Vki_tissue_RE",
            elementRef="Vki_tissue",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vki{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vki_blood_RE",
            metaId="Vki_blood_RE",
            elementRef="Vki_blood",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vurine_RE",
            metaId="Vurine_RE",
            elementRef="Vurine",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"Vurine{PORT_SUFFIX}",
        ),
        # liver
        ReplacedElement(
            sid="Vli_tissue_RE",
            metaId="Vli_tissue_RE",
            elementRef="Vli_tissue",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vli{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vli_blood_RE",
            metaId="Vli_blood_RE",
            elementRef="Vli_blood",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgulumen_RE1",
            metaId="Vgulumen_RE1",
            elementRef="Vgulumen",
            submodelRef=SUBMODEL_SID_DICT["li"],
            portRef=f"Vlumen{PORT_SUFFIX}",
        ),
        # intestine
        ReplacedElement(
            sid="Vstomach_RE",
            metaId="Vstomach_RE",
            elementRef="Vstomach",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vstomach{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgu_blood_RE",
            metaId="Vgu_blood_RE",
            elementRef="Vgu_blood",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vext{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgu_RE",
            metaId="Vgu_RE",
            elementRef="Vgu",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vgu{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vgulumen_RE2",
            metaId="Vgulumen_RE2",
            elementRef="Vgulumen",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vlumen{PORT_SUFFIX}",
        ),
        ReplacedElement(
            sid="Vfeces_RE",
            metaId="Vfeces_RE",
            elementRef="Vfeces",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"Vfeces{PORT_SUFFIX}",
        ),
    ]
)

# -------------------------------------------------------------------------------------------------
# Species
# -------------------------------------------------------------------------------------------------
_m.species = []
for sid, sdict in SUBSTANCES_BODY.items():
    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}"
        else:
            # blood compartment
            sid_ex = f"C{cid}_blood_{sid}"
            cid_ex = f"V{cid}_blood"

        _m.species.append(
            Species(
                sid_ex,
                metaId=f"meta_{sid_ex}",
                initialConcentration=sdict["cinit"],
                compartment=cid_ex,
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} ({cname} blood)",
                hasOnlySubstanceUnits=False,
                port=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

    # tissue species
    for cid, cname in COMPARTMENTS_BODY.items():
        if sdict.get("tissue_distribution") and cid not in [
            "ve",
            "ar",
            "li",
            "ki",
            "gu",
            "po",
            "hv",
        ]:
            sid_ex = f"C{cid}_{sid}"
            cid_ex = f"V{cid}_tissue"

            _m.species.extend(
                [
                    Species(
                        sid_ex,
                        metaId=f"meta_{sid_ex}",
                        initialConcentration=sdict["cinit"],
                        compartment=cid_ex,
                        substanceUnit=sdict["unit"],
                        name=f"{sdict['name']} ({cname})",
                        hasOnlySubstanceUnits=False,
                        port=True,
                        sboTerm=SBO.SIMPLE_CHEMICAL,
                        annotations=annotations.species[sid],
                    )
                ]
            )
            if cid == "re":
                # the slow binding of rapamycin in the tissue
                _m.species.append(
                    Species(
                        f"C{cid}_deep_{sid}",
                        metaId=f"meta_C{cid}_deep_{sid}",
                        initialConcentration=sdict["cinit"],
                        compartment=cid_ex,
                        substanceUnit=sdict["unit"],
                        name=f"{sdict['name']} bound ({cname})",
                        hasOnlySubstanceUnits=False,
                        sboTerm=SBO.SIMPLE_CHEMICAL,
                        annotations=annotations.species[sid],
                    )
                )

    # urine metabolites dummy species
    if sid in ["rx"]:
        _m.species.extend(
            [
                Species(
                    f"Aurine_{sid}",
                    metaId=f"meta_Aurine_{sid}",
                    initialConcentration=0,
                    compartment="Vurine",
                    substanceUnit=sdict["unit"],
                    name=f"{sdict['name']} (urine)",
                    hasOnlySubstanceUnits=True,
                    sboTerm=SBO.SIMPLE_CHEMICAL,
                    annotations=annotations.species[sid],
                ),
            ]
        )

    # feces
    if sid in ["rap", "rx"]:
        _m.species.append(
            Species(
                f"Afeces_{sid}",
                metaId=f"meta_Afeces_{sid}",
                initialConcentration=0,
                compartment="Vfeces",
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} (feces)",
                hasOnlySubstanceUnits=True,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )
    if sid in ["rap", "rx"]:
        _m.species.append(
            Species(
                f"Clumen_{sid}",
                metaId=f"meta_Clumen_{sid}",
                initialConcentration=0,
                compartment="Vgulumen",
                substanceUnit=sdict["unit"],
                name=f"{sdict['name']} (lumen)",
                hasOnlySubstanceUnits=False,
                sboTerm=SBO.SIMPLE_CHEMICAL,
                annotations=annotations.species[sid],
            )
        )

# replace species
replaced_species = {
    "li": ["rap", "rx"],
    "ki": ["rx"],
    "gu": ["rap"],
}

# blood to external species
for tkey, skey_list in replaced_species.items():
    for skey in skey_list:
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"C{tkey}_blood_{skey}_RE",
                    metaId=f"C{tkey}_blood_{skey}_RE",
                    elementRef=f"C{tkey}_blood_{skey}",
                    submodelRef=SUBMODEL_SID_DICT[tkey],
                    portRef=f"{skey}_ext{PORT_SUFFIX}",
                ),
            ]
        )

# urine species replacements
for skey in ["rx"]:
    _m.replaced_elements.append(
        ReplacedElement(
            sid=f"Aurine_{skey}_RE",
            metaId=f"Aurine_{skey}_RE",
            elementRef=f"Aurine_{skey}",
            submodelRef=SUBMODEL_SID_DICT["ki"],
            portRef=f"{skey}_urine{PORT_SUFFIX}",
        )
    )

# feces species replacements
for skey in ["rap", "rx"]:
    _m.replaced_elements.append(
        ReplacedElement(
            sid=f"Afeces_{skey}_RE",
            metaId=f"Afeces_{skey}_RE",
            elementRef=f"Afeces_{skey}",
            submodelRef=SUBMODEL_SID_DICT["gu"],
            portRef=f"{skey}_feces{PORT_SUFFIX}",
        )
    )

# intestine species replacements
for skey in ["rap", "rx"]:
    _m.replaced_elements.extend(
        [
            ReplacedElement(
                sid=f"Clumen_{skey}_RE1",
                metaId=f"Clumen_{skey}_RE1",
                elementRef=f"Clumen_{skey}",
                submodelRef=SUBMODEL_SID_DICT["gu"],
                portRef=f"{skey}_lumen{PORT_SUFFIX}",
            )
        ]
    )

for skey in ["rx"]:
    _m.replaced_elements.extend(
        [
            ReplacedElement(
                sid=f"Clumen_{skey}_RE2",
                metaId=f"Clumen_{skey}_RE2",
                elementRef=f"Clumen_{skey}",
                submodelRef=SUBMODEL_SID_DICT["li"],
                portRef=f"{skey}_lumen{PORT_SUFFIX}",
            ),
        ]
    )

# -------------------------------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------------------------------
_m.parameters.extend(
    [
        # whole body data
        Parameter(
            "BW",
            75,
            U.kg,
            constant=True,
            name="body weight [kg]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "BW_ref",
            75,
            U.kg,
            constant=True,
            name="reference body weight [kg]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=(
                "Body weight at which the exchange flows of the tissue "
                "distribution (`ftissue_<sid>`) are given; they scale with the "
                "body weight as the volumes, the blood flows and the renal "
                "clearance do."
            ),
        ),
        Parameter(
            "HEIGHT",
            170,
            U.cm,
            constant=True,
            name="height [cm]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "HR",
            70,
            U.per_min,
            constant=True,
            name="heart rate [1/min]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "HRrest",
            70,
            U.per_min,
            constant=True,
            name="resting heart rate [1/min]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "BSA",
            0,
            U.m2,
            constant=False,
            name="body surface area [m^2]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "COBW",
            1.548,
            U.ml_per_s_kg,
            constant=True,
            name="cardiac output per bodyweight [ml/s/kg]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cardiac_function",
            1.0,
            U.dimensionless,
            constant=False,
            name="heart function",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            1.0: normal function;
            <1.0: reduced function as reduced cardiac output;
            """,
        ),
        Parameter(
            "CO",
            108.33,
            U.ml_per_s,
            constant=False,
            name="cardiac output [ml/s]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "QC",
            108.33 / 1000 * 60,
            U.l_per_min,
            constant=False,
            name="cardiac output [l/min]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "COHRI",
            150,
            U.ml,
            constant=True,
            name="increase of cardiac output per heartbeat [ml]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue volumes
        Parameter(
            "Fblood",
            0.02,
            U.dimensionless,
            constant=False,
            name="blood fraction of organ volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue volumes
        Parameter(
            "FVgu",
            0.0171,
            U.l_per_kg,
            constant=True,
            name="gut fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVki",
            0.0044,
            U.l_per_kg,
            constant=True,
            name="kidney fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVli",
            0.0210,
            U.l_per_kg,
            constant=True,
            name="liver fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVlu",
            0.0076,
            U.l_per_kg,
            constant=True,
            name="lung fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVre",
            0,
            U.l_per_kg,
            constant=False,
            name="rest of body fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # calculated based on other tissues
        Parameter(
            "FVve",
            0.0514,
            U.l_per_kg,
            constant=True,
            name="venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVar",
            0.0257,
            U.l_per_kg,
            constant=True,
            name="arterial fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVpo",
            0.001,
            U.l_per_kg,
            constant=True,
            name="portal fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FVhv",
            0.001,
            U.l_per_kg,
            constant=True,
            name="hepatic venous fractional tissue volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        # fractional tissue blood flows
        Parameter(
            "FQgu",
            0.18,
            U.dimensionless,
            constant=True,
            name="gut fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
        Assumed FQgu is the sum of FQgu + FQsp + FQpa for correct portal blood flow.
        with FQgu=0.146, FQsp=0.017, FQpa=0.017
        """,
        ),
        Parameter(
            "FQki",
            0.190,
            U.dimensionless,
            constant=True,
            name="kidney fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQh",
            0.215,
            U.dimensionless,
            constant=True,
            name="hepatic (venous side) fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQlu",
            1,
            U.dimensionless,
            constant=True,
            name="lung fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "FQre",
            0,
            U.dimensionless,
            constant=False,
            name="rest of body fractional tissue blood flow",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "f_cirrhosis",
            0,
            U.dimensionless,
            constant=True,
            name="severity of cirrhosis [0, 0.95]",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (BQB.IS, "ncit/C2951"),  # Cirrhosis
                (BQB.IS, "efo/0001422"),  # cirrhosis of liver
            ],
        ),
        Parameter(
            "f_shunts",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of portal venous blood shunted by the liver",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            annotations=[
                (
                    BQB.IS,
                    "hp/HP:0002629",
                ),  # Gastrointestinal arteriovenous malformation
            ],
        ),
        Parameter(
            "f_tissue_loss",
            0,
            U.dimensionless,
            constant=False,
            name="fraction of lost parenchymal liver volume",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
        Parameter(
            "r_cyp3a5",
            **FITS.parameter("r_cyp3a5", unit=U.dimensionless),
            constant=True,
            name="capacity of CYP3A5 relative to CYP3A4",
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes="""
            Metabolism of rapamycin by CYP3A5 of the genotype *1/*1 relative
            to the one by CYP3A4, in the gut wall and the liver. The reference
            genotype is CYP3A5 *3/*3 (no functional CYP3A5, `f_cyp3a5` = 0),
            the genotype of most Caucasians.
            """,
        ),
    ]
)


for sid, sdict in SUBSTANCES_BODY.items():
    if "PODOSE" in sdict:
        # Parameter for intestine model
        _m.parameters.extend(
            [
                Parameter(
                    f"PODOSE_{sid}",
                    sdict["PODOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"oral dose {sid} [mg]",
                ),
            ]
        )
        _m.replaced_elements.extend(
            [
                ReplacedElement(
                    sid=f"PODOSE_{sid}_RE",
                    metaId=f"PODOSE_{sid}_RE",
                    elementRef=f"PODOSE_{sid}",
                    submodelRef=SUBMODEL_SID_DICT["gu"],
                    portRef=f"PODOSE_{sid}{PORT_SUFFIX}",
                ),
            ]
        )

_m.rules = [
    *_m.rules,
    AssignmentRule("f_shunts", "f_cirrhosis", unit=U.dimensionless),
    AssignmentRule("f_tissue_loss", "f_cirrhosis", unit=U.dimensionless),
]

replaced_parameters = {
    "li": ["r_cyp3a5"],
    "gu": ["Mr_rap", "r_cyp3a5"],
    "ki": [],
}

for ckey, pids in replaced_parameters.items():
    for pid in pids:
        _m.replaced_elements.append(
            ReplacedElement(
                sid=f"{pid}_{ckey}_RE",
                metaId=f"{pid}_{ckey}_RE",
                elementRef=f"{pid}",
                submodelRef=SUBMODEL_SID_DICT[ckey],
                portRef=f"{pid}{PORT_SUFFIX}",
            )
        )

# species specific parameters
for sid, sdict in SUBSTANCES_BODY.items():
    if "Mr" in sdict:
        _m.parameters.append(
            Parameter(
                f"Mr_{sid}",
                sdict["Mr"],
                U.g_per_mole,
                constant=True,
                name=f"Molecular weight {sid} [g/mole]",
                sboTerm=SBO.MOLECULAR_MASS,
            )
        )

    # tissue distribution
    if sdict.get("tissue_distribution"):
        _m.parameters.extend(
            [
                Parameter(
                    f"ftissue_{sid}",
                    **FITS.parameter(f"ftissue_{sid}", unit=U.l_per_min),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue distribution {sid}",
                    notes="Exchange flow at the reference body weight `BW_ref`.",
                ),
                Parameter(
                    f"Kp_{sid}",
                    **FITS.parameter(f"Kp_{sid}", unit=U.dimensionless),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"tissue/blood partition coefficient {sid}",
                    notes="""
                    Partition coefficient between the tissue and whole blood;
                    with the exchange flow `ftissue` it sets the volume of
                    distribution (about 12 l/kg of whole blood) and the long
                    terminal half-life of about 60 hr.
                    """,
                ),
                Parameter(
                    f"Kpdeep_{sid}",
                    **FITS.parameter(f"Kpdeep_{sid}", unit=U.dimensionless),
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"bound/free partition coefficient {sid} (tissue)",
                    notes="""
                    Ratio of the bound to the free rapamycin of the tissue in
                    the equilibrium of the slow binding (FKBP12): the volume of
                    distribution of the terminal phase.
                    """,
                ),
                Parameter(
                    f"kdeep_{sid}",
                    **FITS.parameter(f"kdeep_{sid}", unit=U.per_min),
                    constant=True,
                    sboTerm=SBO.KINETIC_CONSTANT,
                    name=f"rate of the slow binding {sid} (tissue)",
                    notes="""
                    The release of the bound rapamycin of the tissue limits the
                    terminal phase: the concentrations fall by a factor of about
                    10 in the first day after a dose (distribution) and with a
                    half-life of about 60 hr after it.
                    """,
                ),
                # blood to plasma ratio, for the plasma concentrations
                Parameter(
                    f"BP_{sid}",
                    sdict["BP"],
                    U.dimensionless,
                    constant=True,
                    name=f"blood to plasma ratio {sid}",
                ),
            ]
        )

    # dosing
    if "IVDOSE" in sdict:
        _m.parameters.extend(
            [
                Parameter(
                    f"IVDOSE_{sid}",
                    sdict["IVDOSE"],
                    U.mg,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"IV bolus dose {sid} [mg]",
                ),
                # iv kinetics after application
                Parameter(
                    f"ti_{sid}",
                    10,
                    U.second,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"injection time {sid} [s]",
                ),
                Parameter(
                    f"Ki_{sid}",
                    0.02,
                    U.per_min,
                    constant=False,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ki [1/min] injection {sid}",
                ),
                # continuous infusion
                Parameter(
                    f"Ri_{sid}",
                    0,
                    U.mg_per_min,
                    constant=True,
                    sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                    name=f"Ri [mg/min] rate of infusion {sid}",
                ),
            ]
        )

# -------------------------------------------------------------------------------------------------
# AssignmentRules
# -------------------------------------------------------------------------------------------------
_m.rules = [
    *_m.rules,
    # Rest body volume
    AssignmentRule(
        "FVre",
        "1.0 l_per_kg - (FVgu + FVki + FVli + FVlu + FVve + FVar)",
        U.l_per_kg,
    ),
    # Rest body perfusion
    AssignmentRule("FQre", "1.0 dimensionless - (FQki + FQh)", U.dimensionless),
    # Body surface area (Haycock1978)
    AssignmentRule(
        "BSA", "0.024265 m2 * power(BW/1 kg, 0.5378) * power(HEIGHT/1 cm, 0.3964)", U.m2
    ),
    # cardiac output (depending on heart rate and bodyweight)
    AssignmentRule(
        "CO",
        "f_cardiac_function*BW*COBW + (HR-HRrest)*COHRI / 60 s_per_min",
        U.ml_per_s,
    ),
    # cardiac output (depending on bodyweight)
    AssignmentRule("QC", "CO/1000 ml_per_l * 60 s_per_min", U.l_per_min),
    # volumes
    AssignmentRule("Vgu", "BW*FVgu", U.liter),
    AssignmentRule("Vki", "BW*FVki", U.liter),
    AssignmentRule("Vli", "BW*FVli", U.liter),
    AssignmentRule("Vlu", "BW*FVlu", U.liter),
    AssignmentRule("Vre", "BW*FVre", U.liter),
    # venous and arterial blood volume (corrected for tissue blood volumes)
    AssignmentRule(
        "Vve",
        "BW*FVve - FVve/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Var",
        "BW*FVar - FVar/(FVar+FVve) * BW * Fblood * (1 l_per_kg - FVve - FVar)",
        U.liter,
    ),
    AssignmentRule(
        "Vpo",
        "BW*FVpo - FVpo/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv))",
        U.liter,
    ),
    AssignmentRule(
        "Vhv",
        "BW*FVhv - FVhv/(FVar+FVve+FVpo+FVhv) "
        "* BW * Fblood * (1 l_per_kg - (FVar+FVve+FVpo+FVhv))",
        U.liter,
    ),
    # blood flows
    AssignmentRule("Qgu", "QC*FQgu", U.l_per_min, name="gut blood flow"),
    AssignmentRule("Qki", "QC*FQki", U.l_per_min, name="kidney blood flow"),
    AssignmentRule(
        "Qh", "QC*FQh", U.l_per_min, name="hepatic (venous side) blood flow"
    ),
    AssignmentRule("Qha", "Qh - Qgu", U.l_per_min, name="hepatic artery blood flow"),
    AssignmentRule("Qlu", "QC*FQlu", U.l_per_min, name="lung blood flow"),
    AssignmentRule("Qre", "QC*FQre", U.l_per_min, name="rest of body blood flow"),
    AssignmentRule("Qpo", "Qgu", U.l_per_min, name="portal blood flow"),
]

# Volumes for explicit tissue models
for cid, cname in COMPARTMENTS_BODY.items():
    if cid not in ["ve", "ar", "po", "hv"]:
        _m.rules.append(
            # blood volume associated with tissue
            AssignmentRule(
                f"V{cid}_blood",
                value=f"V{cid} * Fblood",
                unit=U.liter,
                name=f"blood volume of {cname}",
            ),
        )

        if cid == "li":
            _m.rules.append(
                # Cirrhosis: Adjustment of fractional liver tissue volume
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - f_tissue_loss) * (1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )
        else:
            _m.rules.append(
                AssignmentRule(
                    f"V{cid}_tissue",
                    value=f"V{cid}*(1 dimensionless - Fblood)",
                    unit=U.liter,
                    name=f"tissue volume of {cname}",
                ),
            )

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    # injection
    if "IVDOSE" in sdict:
        _m.rules.extend(
            [
                AssignmentRule(
                    f"Ki_{sid}",
                    f"0.693 dimensionless/ti_{sid} * 60 s_per_min",
                    U.per_min,
                    name="injection rate IV",
                    notes="""
                log(2) ~ 0.693 used for conversion of half-life injection time to
                injection rate.
                """,
                ),
            ]
        )

# --------------------------------------------------------------------------------------------------
# Reactions
# --------------------------------------------------------------------------------------------------
_m.reactions = []

for sid, sdict in SUBSTANCES_BODY.items():
    sname = sdict["name"]

    # --------------------
    # tissue distribution
    # --------------------
    if sdict.get("tissue_distribution"):
        for cid in COMPARTMENTS_BODY:
            if cid not in ["ve", "ar", "li", "gu", "ki", "po", "hv"]:
                _m.reactions.append(
                    Reaction(
                        sid=f"transport_{cid}_{sid}",
                        name=f"transport {sname}",
                        formula=(
                            f"ftissue_{sid} * BW/BW_ref * (C{cid}_blood_{sid}*Kp_{sid} - C{cid}_{sid})",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} <-> C{cid}_{sid}",
                        compartment=f"V{cid}_tissue",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                )
                if cid == "re":
                    _m.reactions.append(
                        Reaction(
                            sid=f"deep_{cid}_{sid}",
                            name=f"slow binding {sname} (tissue)",
                            formula=(
                                f"kdeep_{sid} * V{cid}_tissue * "
                                f"(Kpdeep_{sid} * C{cid}_{sid} - C{cid}_deep_{sid})",
                                U.mmole_per_min,
                            ),
                            equation=f"C{cid}_{sid} <-> C{cid}_deep_{sid}",
                            compartment=f"V{cid}_tissue",
                            sboTerm=SBO.BIOCHEMICAL_REACTION,
                        ),
                    )

    if "IVDOSE" in sdict:
        _m.reactions.extend(
            [
                # --------------------
                # iv dose application
                # --------------------
                Reaction(
                    sid=f"iv_{sid}",
                    name=f"iv {sname}",
                    formula=(f"Ki_{sid}*IVDOSE_{sid}/Mr_{sid}", U.mmole_per_min),
                    equation=f"-> Cve_{sid}",
                    sboTerm=SBO.TRANSPORT_REACTION,
                    annotations=[(BQB.IS, "ncit/C38276")],
                    compartment="Vve",
                ),
            ]
        )

    for cid, cname in COMPARTMENTS_BODY.items():
        if cid in ["ve", "ar", "po", "hv"]:
            continue
        # --------------------
        # ve -> lung -> ar
        # --------------------
        if cid == "lu":
            rid_in = f"Flow_ve_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"

            rid_out = f"Flow_{cid}_ar_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Cve_{sid}", U.mmole_per_min),
                        equation=f"Cve_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Car_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> ve
        # --------------------
        if cid in ["ki", "re"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_ve_{sid}"
            name_out = f"outflow {cname} {sname}"

            # only distribution in blood volume
            _m.reactions.extend(
                [
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> organ -> po
        # --------------------
        if cid in ["gu"]:
            rid_in = f"Flow_ar_{cid}_{sid}"
            name_in = f"inflow {cname} {sname}"
            rid_out = f"Flow_{cid}_po_{sid}"
            name_out = f"outflow {cname} {sname}"

            _m.reactions.extend(
                [
                    # only distribution in blood volume
                    Reaction(
                        sid=rid_in,
                        name=name_in,
                        formula=(f"Q{cid}*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> C{cid}_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=rid_out,
                        name=name_out,
                        formula=(
                            f"Q{cid}*C{cid}_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"C{cid}_blood_{sid} -> Cpo_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )
        # --------------------
        # ar -> li
        # po -> li
        # li -> hv -> ve
        # --------------------
        if cid == "li":
            # liver
            _m.reactions.extend(
                [
                    Reaction(
                        sid=f"Flow_arli_li_{sid}",
                        name=f"arterial inflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qha*Car_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Car_{sid} -> Cli_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # shunted arterial flow
                    Reaction(
                        sid=f"Flow_arli_hv_{sid}",
                        name="flow arterial shunts",
                        formula=(f"f_shunts*Qha*Car_{sid}", U.mmole_per_min),
                        equation=f"Car_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # (unshunted) portal flow
                    Reaction(
                        sid=f"Flow_po_li_{sid}",
                        name=f"outflow po {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*Qpo*Cpo_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cpo_{sid} -> Cli_blood_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    # portal shunts
                    Reaction(
                        sid=f"Flow_po_hv_{sid}",
                        name="flow portal shunts",
                        formula=(f"f_shunts*Qpo*Cpo_{sid}", U.mmole_per_min),
                        equation=f"Cpo_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_li_hv_{sid}",
                        name=f"outflow liver {sname}",
                        formula=(
                            f"(1 dimensionless - f_shunts)*(Qpo+Qha)*Cli_blood_{sid}",
                            U.mmole_per_min,
                        ),
                        equation=f"Cli_blood_{sid} -> Chv_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                    Reaction(
                        sid=f"Flow_hv_ve_{sid}",
                        name=f"outflow hepatic vein {sname}",
                        formula=(f"Qh*Chv_{sid}", U.mmole_per_min),
                        equation=f"Chv_{sid} -> Cve_{sid}",
                        sboTerm=SBO.TRANSPORT_REACTION,
                    ),
                ]
            )

# --------------------------------------------------------------------------------------------------
# RateRules
# --------------------------------------------------------------------------------------------------
_m.rate_rules = []

for sid, sdict in SUBSTANCES_BODY.items():
    if "IVDOSE" in sdict:
        _m.rate_rules.extend(
            [
                # injection of dose
                RateRule(
                    f"IVDOSE_{sid}", f"-iv_{sid}*Mr_{sid} + Ri_{sid}", U.mg_per_min
                ),
            ]
        )

# --------------------------------------------------------------------------------------------------
# Blood and plasma concentrations
# --------------------------------------------------------------------------------------------------
# The circulating species are whole blood concentrations, the measured matrix of
# nearly all studies. Rapamycin is bound in the erythrocytes (FKBP12), its
# plasma concentration is the blood concentration over the blood to plasma
# ratio `BP_rap`. The metabolites are much less bound in the erythrocytes (the
# blood to plasma ratio of the total radioactivity is 2.7 against 142 of
# rapamycin, Leung2006); their plasma concentration is taken as their blood
# concentration.
_m.species.extend(
    [
        Species(
            "Cveblood_rap",
            metaId="meta_Cveblood_rap",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name="rapamycin (venous blood)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["rap"],
        ),
        Species(
            "Cveplasma_rap",
            metaId="meta_Cveplasma_rap",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name="rapamycin (venous plasma)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["rap"],
        ),
        Species(
            "Cveblood_raptot",
            metaId="meta_Cveblood_raptot",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name="rapamycin and metabolites (venous blood)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["raptot"],
        ),
        Species(
            "Cveplasma_raptot",
            metaId="meta_Cveplasma_raptot",
            initialConcentration=0,
            compartment="Vve",
            substanceUnit=U.mmole,
            name="rapamycin and metabolites (venous plasma)",
            hasOnlySubstanceUnits=False,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["raptot"],
        ),
        Species(
            "Afeces_raptot",
            metaId="meta_Afeces_raptot",
            initialConcentration=0,
            compartment="Vfeces",
            substanceUnit=U.mmole,
            name="rapamycin and metabolites (feces)",
            hasOnlySubstanceUnits=True,
            sboTerm=SBO.SIMPLE_CHEMICAL,
            annotations=annotations.species["raptot"],
        ),
    ]
)

_m.rules.extend(
    [
        AssignmentRule(
            "Cveblood_rap",
            "Cve_rap",
            unit=U.mM,
            name="rapamycin blood concentration",
        ),
        AssignmentRule(
            "Cveplasma_rap",
            "Cve_rap / BP_rap",
            unit=U.mM,
            name="rapamycin plasma concentration",
            notes="The blood concentration over the blood to plasma ratio (BP).",
        ),
        AssignmentRule(
            "Cveblood_raptot",
            "Cve_rap + Cve_rx",
            unit=U.mM,
            name="rapamycin and metabolites blood concentration",
            notes="For the comparison with the total radioactivity of blood.",
        ),
        AssignmentRule(
            "Cveplasma_raptot",
            "Cve_rap / BP_rap + Cve_rx",
            unit=U.mM,
            name="rapamycin and metabolites plasma concentration",
            notes="For the comparison with the total radioactivity of plasma.",
        ),
        AssignmentRule(
            "Afeces_raptot",
            "Afeces_rap + Afeces_rx",
            unit=U.mmole,
            name="rapamycin and metabolites in the feces",
        ),
    ]
)

model_body: Model = _m


if __name__ == "__main__":
    from pkdb_models.models.rapamycin import MODEL_BASE_PATH

    result = create_model(
        filepath=MODEL_BASE_PATH / f"{model_body.sid}.xml",
        model=model_body,
        sbml_level=3,
        sbml_version=2,
    )
