"""Sensitivity analysis of the rapamycin model.

The outputs are the pharmacokinetic parameters of rapamycin in whole blood
after a single oral dose (the non-compartmental analysis of `pkpdutils`). The
parameters are the parameters of the model with bounds of +-15 % around their
value. The analyses are run as a script, the local analysis by default:

    python -m pkdb_models.models.rapamycin.sensitivity.sensitivity_analysis
    python -m pkdb_models.models.rapamycin.sensitivity.sensitivity_analysis local sampling

The analyses run on `pkpdsim.sensitivity.SensitivityRunner`, which writes the
table of the parameters into `results/sensitivity` of the model and the tables
and figures of an analysis into `results/sensitivity/<analysis>`.
An analysis does not vary a parameter which one of its groups sets
(`pkpdsim.sensitivity.free_parameters`): `sampling`, over the groups, leaves out
`f_cirrhosis`. The environment variable `PKPDSIM_SENSITIVITY_N` replaces the
numbers of samples of the global analyses for a quick run, see
`pkpdsim.sensitivity.sample_size`.

The model without rapamycin is stationary, so the simulation of a group starts
with the dose. The dose is followed for `END_HR` hours, long enough for the
terminal phase of rapamycin (a half-life of about 60 hr).
"""

from collections.abc import Mapping

import numpy as np
from pkpdsim.sensitivity import (
    FastAnalysis,
    LocalAnalysis,
    ParameterType,
    SamplingAnalysis,
    SensitivityGroup,
    SensitivityOutput,
    SensitivityParameter,
    SensitivityRunner,
    SobolAnalysis,
    parameters_from_sbml,
    set_bounds,
    set_relative_bounds,
)
from pkpdutils import NCAOptions
from sbmlsim import Q
from sbmlsim.simulation import PK, Formula, Simulation

from pkdb_models.models.rapamycin import MODEL, MODEL_PATH
from pkdb_models.models.rapamycin.experiments.base_experiment import (
    RapamycinSimulationExperiment,
)
from pkdb_models.models.rapamycin.fitting.parameters import (
    parameters_pk as fit_parameters,
)

#: oral dose of rapamycin of the analysis [mg]
DOSE_RAP = 5
#: the hours after the dose
END_HR = 336
#: the limit of quantification of the non-compartmental analysis [mM], about
#: 1e-4 of the maximal blood concentration after 5 mg (3e-8 mM), so that the
#: noise of the integrator is not the terminal phase
LLOQ = 3e-12
#: fraction of the bounds of a parameter relative to its value
BOUNDS_FRACTION = 0.15

#: outputs of the non-compartmental analysis of rapamycin -> unit of the output
PK_OUTPUTS: dict[str, str] = {
    "auc_inf_obs": "mM*min",
    "cmax": "mM",
    "tmax": "min",
    "thalf": "min",
    "vz_f": "l",
    "cl_f": "ml/min",
}
#: the outputs whose unit has the amount of the dose: the dose is a mass
#: (`PODOSE_rap` in mg), so they are divided by the molecular weight (its value
#: at the start, a value per simulation like the parameter of the analysis)
PK_PER_DOSE_MASS: frozenset[str] = frozenset({"vz_f", "cl_f"})

#: the outputs with their ids and names
OUTPUTS: list[SensitivityOutput] = [
    SensitivityOutput("rap_auc_inf_obs", "rapamycin AUC∞"),
    SensitivityOutput("rap_cmax", "rapamycin Cmax"),
    SensitivityOutput("rap_tmax", "rapamycin tmax"),
    SensitivityOutput("rap_thalf", "rapamycin Half-life"),
    SensitivityOutput("rap_vz_f", "rapamycin Vz/F"),
    SensitivityOutput("rap_cl_f", "rapamycin CL/F"),
]


def _pk_observables(sid: str, outputs: Mapping[str, str]) -> list:
    """Get the non-compartmental analysis of a substance in blood and its outputs."""
    pk_id = f"pk_{sid}"
    return [
        PK(
            pk_id,
            f"[Cveblood_{sid}]",
            dose=f"PODOSE_{sid}",
            route="oral",
            options=NCAOptions(lloq=LLOQ),
            parameters=list(outputs),
        ),
        *(
            Formula(
                f"{sid}_{key}",
                f"{pk_id}.{key} / at(Mr_{sid}, 0)"
                if key in PK_PER_DOSE_MASS
                else f"{pk_id}.{key}",
                unit=unit,
            )
            for key, unit in outputs.items()
        ),
    ]


#: the observables of a simulation, the outputs with their units; the id of
#: the analysis (`pk_rap`) is no species of the model
OBSERVABLES = _pk_observables("rap", PK_OUTPUTS)

#: the output times [hr]: every 0.02 hr in the first day after the dose
#: (absorption and maximum), every 0.5 hr after it
TIMES_HR: tuple[float, ...] = tuple(
    np.unique(
        np.concatenate(
            [
                np.linspace(0, 24, 1201),
                np.linspace(24, END_HR, 2 * (END_HR - 24) + 1),
            ]
        )
    ).tolist()
)


def sensitivity_simulation(changes: Mapping[str, float]) -> Simulation:
    """Get the simulation of a group: a single oral dose of rapamycin at time 0."""
    return Simulation(
        start=0,
        end=END_HR,
        times=TIMES_HR,
        time_unit="hr",
        preinit_changes={**changes, "PODOSE_rap": Q(DOSE_RAP, "mg")},
    )


#: the groups the analyses are run for: the control and the stages of
#: cirrhosis of the simulation experiments
SENSITIVITY_GROUPS: list[SensitivityGroup] = [
    SensitivityGroup("control", "Control", sensitivity_simulation({}), "dimgrey"),
    *(
        SensitivityGroup(
            uid,
            name,
            sensitivity_simulation(
                {"f_cirrhosis": RapamycinSimulationExperiment.cirrhosis_map[name]}
            ),
            RapamycinSimulationExperiment.cirrhosis_colors[name],
        )
        for uid, name in [
            ("mildCirrhosis", "Mild cirrhosis"),
            ("moderateCirrhosis", "Moderate cirrhosis"),
            ("severeCirrhosis", "Severe cirrhosis"),
        ]
    ),
]


def sensitivity_parameters() -> list[SensitivityParameter]:
    """Get the parameters of the analyses with their bounds and references."""
    parameters = parameters_from_sbml(
        MODEL_PATH,
        exclude_ids={
            # molecular weights
            "Mr_rap",
            # reference values
            "BW_ref",
            # unchangeable values
            "FQlu",
            "FVhv",
            "FVpo",
            # dosing parameters
            "PODOSE_rap",
            "IVDOSE_rap",
            "ti_rap",
            # unused volumes
            "Vurine",
            "Vfeces",
            "Vstomach",
            "GU__Vchain",
            # the fast equilibria of the liver
            "LI__RAPIM_k",
            "LI__RXEX_k",
        },
        exclude_na=True,
        exclude_zero=True,
    )

    # the fitted parameters, their bounds are set relative to the value below
    set_bounds(
        parameters,
        [(fp.pid, np.nan, np.nan, ParameterType.FIT) for fp in fit_parameters],
    )
    # the scaling factors of the organ functions
    set_bounds(
        parameters,
        [
            (uid, 1 - BOUNDS_FRACTION, 1 + BOUNDS_FRACTION, ParameterType.SCALING)
            for uid in ["f_cardiac_function", "KI__f_renal_function"]
        ],
    )

    # the references of the values from the literature
    references = {
        "BW": r"\cite{Ogden2004, Jones2013, Thompson2009, Brown1997}",
        "FQgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FQh": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FQki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVar": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVgu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVki": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVli": r"\cite{Jones2013, Wynne1989, Thompson2009, Brown1997}",
        "FVlu": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "FVve": r"\cite{Jones2013, Thompson2009, Brown1997}",
        "COBW": r"\cite{Cattermole2017, Patel2021, Collis2001}",
        "BP_rap": r"\cite{Leung2006}",
    }
    by_uid = {p.uid: p for p in parameters}
    for uid, reference in references.items():
        p = by_uid[uid]
        p.reference = reference
        if p.type == ParameterType.NA:
            p.type = ParameterType.DATA

    # every parameter without bounds varies around its value
    set_relative_bounds(parameters, BOUNDS_FRACTION)
    return parameters


#: the analyses of the script, with their groups and numbers of samples
SENSITIVITY = SensitivityRunner(
    model=MODEL,
    parameters=sensitivity_parameters,
    outputs=OUTPUTS,
    observables=OBSERVABLES,
    analyses={
        "local": LocalAnalysis(groups=SENSITIVITY_GROUPS[:1], delta=0.01),
        "sampling": SamplingAnalysis(groups=SENSITIVITY_GROUPS, n=1000),
        "fast": FastAnalysis(groups=SENSITIVITY_GROUPS[:1], n=1000),
        "sobol": SobolAnalysis(groups=SENSITIVITY_GROUPS[:1], n=4096),
    },
    default=("local",),
)


if __name__ == "__main__":
    SENSITIVITY.main()
