"""Parameter fit problems for rapamycin.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the training data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The validation data has three parts, each a different question: the
interactions (`INTERACTION_VALIDATION`, the arms with a coadministered drug
whose interaction the model describes with factors which are not fitted), the
doses (`DOSE_VALIDATION`, the repeated doses of the studies of a steady state,
which the single doses of the training predict) and the limitations
(`LIMITATION_VALIDATION`, data of a population, a formulation or a process the
model does not describe). The summaries of data which is fitted one by one are
excluded (`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import mapping_metadata
from pkpdsim.metadata import Coadministration
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.rapamycin import DATA_PATHS, RAPAMYCIN_PATH
from pkdb_models.models.rapamycin.experiments.studies import (
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
)

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share, so a named mapping is named with
#: its study
type NamedMappings = dict[str, set[str]]

EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
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
]

#: coadministered drugs without an interaction with rapamycin: their arms are
#: training data like the arms without a coadministration
NO_INTERACTION: frozenset[Coadministration] = frozenset(
    {
        Coadministration.NONE,
        Coadministration.PREDNISONE,
        Coadministration.MYCOPHENOLATE_MOFETIL,
        Coadministration.FAMOTIDINE,
    }
)

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {}

#: the validation of the interactions: the arms with a coadministered drug
#: which interacts with rapamycin, simulated with factors of CYP3A4, CYP3A5
#: and P-glycoprotein which are not fitted. `filter_coadministration` lets the
#: validation data through for this reason, every other interaction stays out
#: of the selection altogether
INTERACTION_VALIDATION: NamedMappings = {
    # diltiazem, gut CYP3A and P-glycoprotein
    "Bottiger2001": {"fm_rap_RAP10_DIL120"},
    # cyclosporine, gut CYP3A and P-glycoprotein, steady state of renal
    # transplant patients
    "Kelly1999": {
        "fm_rap_RAP_L__CSA",
        "fm_rap_RAP_T2__CSA",
        "fm_rap_RAP_T4__CSA",
        "fm_rap_RAP_T8__CSA",
    },
    # tacrolimus, no interaction in the data (AUC ratios 0.97 and 1.07)
    "Tortorici2013": {"fm_rap_TAC10", "fm_rap_TAC10_4H"},
    # rifampin, induction of gut CYP3A4 and P-glycoprotein
    "Tortorici2014": {"fm_rap_RAP20_RIF600"},
    # ritonavir, inhibition of CYP3A4 and CYP3A5
    "Zha2022": {"fm_rap_RAP05_3D"},
    # cyclosporine together with rapamycin and 4 hr before
    "Zimmerman2003": {"fm_rap_RAP10_CSA300", "fm_rap_RAP10_CSA300B"},
    # cyclosporine of renal transplant patients, repeated doses
    "Zimmerman1997": {
        "fm_rap_RAP05",
        "fm_rap_RAP1_5",
        "fm_rap_RAP2_5",
        "fm_rap_RAP3_5",
        "fm_rap_RAP6_5",
    },
}

#: the validation of the doses: the repeated doses of the studies of a steady
#: state, which the single doses of the training predict
DOSE_VALIDATION: NamedMappings = {}

#: the validation of the limitations of the model: valid data of a population,
#: a formulation or a process the model does not describe, which would distort
#: the fitted parameters of the processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    # a nano-amorphous powder, not the formulation of Rapamune: its exposure
    # per dose is 2-6 times the one of Rapamune
    "BasaDenes2019": {
        "fm_rap_RAP05",
        "fm_rap_RAP2",
        "fm_rap_RAP10",
        "fm_rap_RAP40",
        "fm_rap_RAP40_FD",
    },
    # the total radioactivity: the blood to plasma ratio of the metabolites is
    # an assumption of the model
    "Leung2006": {"fm_raptot_Fig1_blood", "fm_raptot_Fig1_plasma"},
    # renal transplant patients on their comedication, the doses of the
    # genotype groups are estimated
    "Renders2007": {"fm_rap_CYP3A5_1_3", "fm_rap_CYP3A5_3_3"},
    # obese patients after a gastric bypass, which the model does not have
    "Rogers2008": {"fm_rap_P5", "fm_rap_P6", "fm_rap_P7", "fm_rap_P8"},
    # renal transplant patients with mixed CYP3A5 genotypes, the doses are
    # derived from the dose normalised values
    "Wang2014": {"fm_rap_RAP", "fm_rap_RAP_CSA"},
    # severe cirrhosis (Child-Pugh C), the shared cirrhosis model
    "Zimmerman2008": {"fm_rap_RAP15_GRADE_C"},
}

#: the summaries of data which the fits use one by one, excluded from every
#: fit: a mean of subjects would count the subjects twice
AGGREGATE_MAPPINGS: NamedMappings = {
    # the 31 subjects of the CYP3A5 groups split by CYP3A4: the *1G activity
    # of the model is derived from this split
    "Zhang2017": {
        "fm_rap_RAP5_CYP3A4_1_1",
        "fm_rap_RAP5_CYP3A4_1_1G",
        "fm_rap_RAP5_CYP3A4_1G_1G",
    },
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(
            mappings.get(study, set())
            for mappings in [
                INTERACTION_VALIDATION,
                DOSE_VALIDATION,
                LIMITATION_VALIDATION,
            ]
        )
    )
    for study in INTERACTION_VALIDATION | DOSE_VALIDATION | LIMITATION_VALIDATION
}


def study(fit_mapping: FitMapping) -> str:
    """Get the study of a fit mapping, the id of its experiment."""
    return fit_mapping.reference.experiment.sid


def named(
    mappings: NamedMappings, fit_mapping_key: str, fit_mapping: FitMapping
) -> bool:
    """Check whether a fit mapping is one of the named mappings of its study."""
    return fit_mapping_key in mappings.get(study(fit_mapping), set())


def check_named(mappings: NamedMappings, fit_mappings: FitMappings) -> None:
    """Check that the named mappings are fit mappings of their studies.

    Raises:
        ValueError: for a study which is no experiment, or a key which is no fit
            mapping of its study.
    """
    for sid, keys in mappings.items():
        experiment = fit_mappings.experiments.get(sid)
        if experiment is None:
            raise ValueError(f"'{sid}' is no experiment of the fits.")
        if unknown := sorted(keys - set(experiment._fit_mappings)):
            raise ValueError(f"'{unknown}' are no fit mappings of '{sid}'.")


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No summaries of data which is fitted one by one, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


def filter_coadministration(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No interacting drug, besides the validation data."""
    if named(VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping):
        return True
    return mapping_metadata(fit_mapping).coadministration in NO_INTERACTION


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the training data, the mappings which fail them are
    excluded, and so are the summaries of data which is fitted one by one
    (`AGGREGATE_MAPPINGS`) and the arms with an interacting drug
    (`filter_coadministration`). The outliers are `OUTLIER_MAPPINGS`, the
    validation data is `VALIDATION_MAPPINGS`. Neither kind enters the cost;
    both are simulated and evaluated with the training data, so a report says
    how the fitted parameters describe data they were not fitted on.

    Args:
        filters: filters of the training data.
        validation: the keys of the validation data by study.

    Returns:
        The fit mapping collections of all kinds by experiment id.

    Raises:
        ValueError: if an outlier is a mapping of no study or of several.
    """
    fit_mappings = FitMappings(
        experiment_classes=EXPERIMENT_CLASSES,
        base_path=RAPAMYCIN_PATH,
        data_path=DATA_PATHS,
    )
    for mappings in [OUTLIER_MAPPINGS, AGGREGATE_MAPPINGS, validation]:
        check_named(mappings, fit_mappings)
    # sbmlsim selects the outliers by their keys alone, which must name one
    # mapping of all studies
    outliers = {key for keys in OUTLIER_MAPPINGS.values() for key in keys}
    for key in outliers:
        studies = [
            sid
            for sid, experiment in fit_mappings.experiments.items()
            if key in experiment._fit_mappings
        ]
        if len(studies) != 1:
            raise ValueError(f"The outlier '{key}' is a mapping of {studies}.")
    selection = fit_mappings.select(
        filters=[
            *(filters if isinstance(filters, list) else [filters]),
            filter_aggregates,
            filter_coadministration,
        ],
        outliers=outliers,
        validation=lambda key, fit_mapping: named(validation, key, fit_mapping),
        print_info=False,
    )
    # the counts per experiment and kind, not the single mappings
    selection.print(detail=False)
    # the validation mappings by the question they answer
    selected = {
        (sid, key)
        for sid, kinds in selection.kinds.items()
        for key, kind in kinds.items()
        if kind is MappingKind.VALIDATION
    }

    def pairs(mappings: NamedMappings) -> set[tuple[str, str]]:
        return {(sid, key) for sid, keys in mappings.items() for key in keys}

    console.print(
        "validation: "
        + ", ".join(
            f"{len(selected & pairs(mappings))} {question}"
            for question, mappings in [
                ("interactions", INTERACTION_VALIDATION),
                ("doses", DOSE_VALIDATION),
                ("limitations", LIMITATION_VALIDATION),
            ]
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """Pharmacokinetic data."""
    return mapping_collections([])


if __name__ == "__main__":
    console.rule(title=f_collections_pk.__name__, align="left", style="white")
    f_collections_pk()
