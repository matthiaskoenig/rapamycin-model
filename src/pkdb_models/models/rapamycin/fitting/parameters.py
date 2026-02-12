"""FitParameters for rapamycin fitting."""

from sbmlsim.fit import FitParameter


parameters_pharmacokinetics = [
    # tissue distribution
    FitParameter(
        pid="ftissue_rap",
        lower_bound=0.01,
        start_value=0.1,
        upper_bound=10,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_rap",
        lower_bound=10,
        start_value=1000,
        upper_bound=10000,
        unit="dimensionless",
    ),


    # absorption
    FitParameter(
        pid="GU__RAPIM_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="GU__RAP2RX_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="GU__RXPG_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="1/min",
    ),

    # excretion feces
    FitParameter(
        pid="GU__RXEXC_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="1/min",
    ),

    # kidney removal
    FitParameter(
        pid="KI__RXEX_k",
        lower_bound=1E-1,
        start_value=1.0,
        upper_bound=10,
        unit="1/min",
    ),

    # hepatic metabolism
    FitParameter(
        pid="LI__RAP2RX_Vmax",
        lower_bound=1E-4,
        start_value=0.02,
        upper_bound=1,
        unit="mmole/min/l",
    ),

    # biliary excretion
    FitParameter(
        pid="LI__RXBEX_k",
        lower_bound=1E-7,
        start_value=0.0001,
        upper_bound=1E-3,
        unit="1/min",
    ),
]


parameters_all = parameters_pharmacokinetics
