"""Annotations of the rapamycin models: the model and its species.

The compartments are the shared ones of `pkpdsim.annotations`.
"""

from pkpdsim.annotations import compartments
from sbmlutils.metadata import (
    BQB,
)

__all__ = ["compartments", "model", "species"]

model = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human
    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model
    # biological process explained by model (GO/NCIT)
    (BQB.IS, "NCIT:C1212"),  # rapamycin/Sirolimus
    (BQB.HAS_PROPERTY, "NCIT:C15261"),  # immunosuppressive therapy
    (BQB.HAS_PROPERTY, "NCIT:C15299"),  # pharmacokinetics
    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "NCIT:C15342"),  # transplantation
    # model reference
    # (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.17091694"),
    # model
    # (BQM.IS, "biomodels.db/MODEL2510140001"),
]

species = {
    "rap": [  # rapamycin
        (BQB.IS, "pubchem.compound/5284616"),
        (BQB.IS, "inchikey/QFJCIRLUMZQUOT-HPLJOQBZSA-N"),
        (BQB.IS, "chebi/CHEBI:9168"),
        (BQB.IS, "ncit/C1212"),
    ],
    "rx": [  # rapamycin metabolite
    ],
    "raptot": [  # rapamycin total metabolites
    ],
}
