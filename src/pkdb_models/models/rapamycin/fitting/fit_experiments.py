"""Subsets of data for rapamycin."""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.rapamycin import RAPAMYCIN_PATH, DATA_PATHS
from pkdb_models.models.rapamycin.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, RapamycinMappingMetaData, Coadministration, Genotype
)
from pkdb_models.models.rapamycin.experiments.studies import *

logger = get_logger(__name__)

# --- Fit experiments ---
f_fitexp_kwargs = dict(
    experiment_classes  = [
        BasaDenes2019,
        Bottiger2001,
        Brattstrom2000,
        Kelly1999,
        KorthBradley2012,
        Leelahavanichkul2005,
        Leung2006,
        Renders2007,
        Rogers2008,
        Tortorici2013,
        Tortorici2014,
        Wang2014,
        Zha2022,
        Zhang2017,
        Zimmerman1997,
        Zimmerman1999,
        Zimmerman2003,
        Zimmerman2008,
    ],
    base_path=RAPAMYCIN_PATH,
    data_path=DATA_PATHS,
)

# --- Filters ---
def filter_baseline(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Healthy control model with wildtype genotype."""

    metadata: RapamycinMappingMetaData = fit_mapping.metadata

    # filter co-administration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, ...)
    if metadata.health not in {Health.HEALTHY}:
        return False

    if metadata.dosing not in {Dosing.SINGLE}:
        return False

    # filter genotype
    if metadata.genotype not in {Genotype.CYP3A4_1_1, Genotype.CYP3A5_1_1, Genotype.NR}:
        return False

    # remove outliers
    if metadata.outlier:
        return False

    return True


def f_fitexp_all():
    """All data."""
    return f_fitexp(metadata_filters=filter_empty, **f_fitexp_kwargs)


def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(metadata_filters=[filter_baseline], **f_fitexp_kwargs)


if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        f_fitexp_all,
        f_fitexp_control,

    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
