"""The fit parameters of the rapamycin model.

The pharmacokinetics are fitted on whole blood concentrations after oral doses.
No study has intravenous data, so the clearance and the bioavailability cannot
be separated: the fraction absorbed (`GU__F_rap_abs`) and the gut wall
(`GU__RAP2RX_k`, `GU__RAPPG_k`, `GU__RAPBASO_k`) are fixed to an oral
bioavailability of about 14 % (Rapamune label), the fit adjusts the hepatic
clearance, the tissue distribution, the dissolution, the effects of food and
the capacity of CYP3A5.
"""

from sbmlsim.fit import FitParameter

#: the pharmacokinetics
parameters_pk: list[FitParameter] = [
    # rate of dissolution and gastric emptying, tmax about 1 hr
    FitParameter(
        pid="GU__Ka_dis_rap",
        start_value=2.0,
        lower_bound=0.1,
        upper_bound=20,
        unit="1/hr",
    ),
    # hepatic clearance of whole blood (about 3 l/hr in total)
    FitParameter(
        pid="LI__RAP2RX_k",
        start_value=0.03,
        lower_bound=1e-3,
        upper_bound=1,
        unit="1/min",
    ),
    # tissue distribution: the fast distribution of the first day
    FitParameter(
        pid="Kp_rap",
        start_value=3.0,
        lower_bound=0.3,
        upper_bound=30,
        unit="dimensionless",
    ),
    FitParameter(
        pid="ftissue_rap",
        start_value=0.5,
        lower_bound=1e-3,
        upper_bound=10,
        unit="l/min",
    ),
    # the slow binding in the tissue: the terminal phase after the
    # distribution of the first day
    FitParameter(
        pid="Kpdeep_rap",
        start_value=5.0,
        lower_bound=0.1,
        upper_bound=100,
        unit="dimensionless",
    ),
    FitParameter(
        pid="kdeep_rap",
        start_value=0.002,
        lower_bound=1e-5,
        upper_bound=0.1,
        unit="1/min",
    ),
    # food: a high fat meal increases the AUC by 35 % and delays tmax
    # (Zimmerman1999)
    FitParameter(
        pid="GU__FED_f_absorption",
        start_value=1.35,
        lower_bound=1.0,
        upper_bound=1.9,
        unit="dimensionless",
    ),
    FitParameter(
        pid="GU__FED_f_dissolution",
        start_value=0.2,
        lower_bound=0.01,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    # CYP3A5 *1/*1 relative to CYP3A4 (Zhang2017, Renders2007)
    FitParameter(
        pid="r_cyp3a5",
        start_value=0.5,
        lower_bound=0.01,
        upper_bound=5.0,
        unit="dimensionless",
    ),
]
