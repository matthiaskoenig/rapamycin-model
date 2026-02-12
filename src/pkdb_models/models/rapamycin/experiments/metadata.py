from dataclasses import dataclass
from enum import Enum

from sbmlsim.fit.objects import MappingMetaData


class Tissue(str, Enum):
    BLOOD = "blood"
    PLASMA = "plasma"
    SERUM = "serum"
    URINE = "urine"
    FECES = "feces"


class Route(str, Enum):
    PO = "po"
    IV = "iv"
    SC = "sc"
    NR = "not reported"


class Dosing(str, Enum):
    SINGLE = "single"
    MULTIPLE = "multiple"
    CONSTANT_INFUSION = "infusion"


class ApplicationForm(str, Enum):
    SUBCUTANEOUS = "subcutaneous"
    TABLET = "tablet"
    SOLUTION = "solution"
    CAPSULE = "capsule"
    MIXED = "mixed"  # mix of forms, e.g. po and iv
    NR = "not reported"


class Health(str, Enum):
    HEALTHY = "healthy"
    HYPERTENSION = "hypertension"
    CIRRHOSIS = "cirrhosis"
    RENAL_IMPAIRMENT = "renal impairment"
    HEPATIC_IMPAIRMENT = "hepatic impairment"
    CHF = "congestive heart failure"
    OBESE = "obese"
    RENAL_TRANSPLANT = "renal transplant"


class Fasting(str, Enum):
    NR = "not reported"
    FASTED = "fasted"
    FED = "fed"


class Coadministration(str, Enum):
    NONE = "none"
    DILTIAZEM = "diltiazem"
    FUJIMYCIN = "fujimycin"  # tacrolimus
    RIFAMPIN = "rifampin"
    CYCLOSPORINE = "cycloposporine"
    RITONAVIR = "ritonavir"  # 3D regimen
    PREDNISONE = "prednisone"
    MPM = "MPM"  # mycophenolic motefil
    FUJIMYCIN_MPM = "fujimycin and mpm"
    FAMOTIDINE = "famotidine"


class Genotype(str, Enum):
    NR = "not reported"

    CYP3A4_1_1 = "CYP3A4 *1/*1",
    CYP3A4_1_1G = "CYP3A4 *1/*1G",
    CYP3A4_1G_1G = "CYP3A4 *1G/*1G",
    CYP3A5_1_1 = "CYP3A5 *1/*1",
    CYP3A5_1_3 = "CYP3A5 *1/*3",
    CYP3A5_3_3 = "CYP3A5 *3/*3",


@dataclass
class RapamycinMappingMetaData(MappingMetaData):
    """Metadata for fitting experiment."""
    tissue: Tissue
    route: Route
    application_form: ApplicationForm
    dosing: Dosing
    health: Health
    fasting: Fasting
    coadministration: Coadministration
    genotype: Genotype = Genotype.NR
    outlier: bool = False

    def to_dict(self):
        return {
            "tissue": self.tissue.name,
            "route": self.route.name,
            "application_form": self.application_form.name,
            "dosing": self.dosing.name,
            "health": self.health.name,
            "fasting": self.fasting.name,
            "coadministration": self.coadministration.name,
            "genotype": self.genotype.name,
            "outlier": self.outlier,
        }
