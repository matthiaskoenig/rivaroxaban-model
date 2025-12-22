"""FitParameters."""

from sbmlsim.fit import FitParameter


parameters_pharmacokinetics = [
    # dissolution rate
    FitParameter(
        pid="GU__Ka_dis_riv",
        lower_bound=1E-3,
        start_value=2,
        upper_bound=100,
        unit="1/hr",
    ),
    # absorption
    FitParameter(
        pid="GU__RIVABS_k",
        lower_bound=0.001,
        start_value=0.1,
        upper_bound=10,
        unit="1/min",
    ),
    # liver metabolism
    FitParameter(
        pid="LI__RIV2RX_Vmax",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=100,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__RXEX_BI_k",
        lower_bound=1E-6,
        start_value=0.0001,
        upper_bound=1E-4,
        unit="1/min",
    ),
    # kidney excretion
    FitParameter(
        pid="KI__RIVEX_k",
        lower_bound=1E-4,
        start_value=1,
        upper_bound=10,
        unit="1/min",
    ),
    # FitParameter(
    #     pid="KI__RXEX_k",
    #     lower_bound=1E-4,
    #     start_value=1,
    #     upper_bound=10,
    #     unit="1/min",
    # ),
]

parameters_pharmacodynamics = [
    FitParameter(
        pid="Emax_PT",
        lower_bound=0.1,
        start_value=1,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(
        pid="EC50_riv_PT",
        lower_bound=1E-7,
        start_value=0.00034,
        upper_bound=1E-2,
        unit="mM",
    ),
    FitParameter(
        pid="Emax_aPTT",
        lower_bound=0.1,
        start_value=1,
        upper_bound=10,
        unit="dimensionless",
    ),
    FitParameter(
        pid="EC50_riv_aPTT",
        lower_bound=1E-7,
        start_value=0.00034,
        upper_bound=1E-2,
        unit="mM",
    ),
    FitParameter(
        pid="Emax_Xa",
        lower_bound=0.5,
        start_value=0.7,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="EC50_riv_Xa",
        lower_bound=1E-7,
        start_value=0.00034,
        upper_bound=1E-2,
        unit="mM",
    ),
]


parameters_all = parameters_pharmacokinetics + parameters_pharmacodynamics