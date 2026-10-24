"""The fit parameters of the rivaroxaban model, with their bounds.

The pharmacokinetics (`parameters_pk`) are fitted first, on the plasma
concentrations and the excretion of rivaroxaban after oral tablets and
solutions: the solutions separate the distribution and the clearance from the
dissolution, the dose range of the tablets (1.25-80 mg) and the fed state
identify the capacity of the dissolution and the effects of food, the urinary
excretion separates the renal from the hepatic clearance and the mass balance
of [14C]rivaroxaban (Weinz2009) the excretion of the metabolites. Fixed with
documented values: the fraction absorbed of the dissolved drug
(`GU__F_riv_abs`, Weinz2009), the unbound fraction (`KI__fu_riv`,
Kubitza2013), the undissolved amount of half-maximal dissolution
(`GU__Km_dis_riv`) and the renal clearance of the metabolites (`KI__RXEX_f`,
not identifiable: the metabolites are cleared fast, their profile fixes only
the split between bile and urine). The pharmacodynamics (`parameters_pd`) are
fitted on the adopted pharmacokinetics.
"""

from sbmlsim.fit import FitParameter

#: the pharmacokinetics
parameters_pk: list[FitParameter] = [
    # absorption of the dissolved drug, the tmax of an oral solution 0.5 hr;
    # it includes the gastric emptying of the solution (half-life 10-20 min),
    # so it is not faster than a half-life of 7 min
    FitParameter(
        pid="GU__RIVABS_k",
        lower_bound=0.005,
        start_value=0.03,
        upper_bound=0.1,
        unit="1/min",
    ),
    # gastric emptying of the tablet, three stages: the tmax of the tablets
    # (2.5-3 hr fasted) is later than the one of the solution (0.5 hr) at
    # doses which dissolve completely (Kubitza2005a)
    FitParameter(
        pid="GU__Ka_emptying_riv",
        lower_bound=0.5,
        start_value=3.0,
        upper_bound=12.0,
        unit="1/hr",
    ),
    # capacity of the dissolution: the Cmax of 10-80 mg tablets
    FitParameter(
        pid="GU__Vmax_dis_riv",
        lower_bound=0.5,
        start_value=7.0,
        upper_bound=50,
        unit="mg/hr",
    ),
    # transit of the undissolved drug: the bioavailability falls with the
    # dose; bounded by the small intestinal transit of 2-7 hr, a faster loss
    # lowers the bioavailability of small doses below the 80-100 % of 10 mg
    # (Mueck2014)
    FitParameter(
        pid="GU__ktransit_riv",
        lower_bound=0.15,
        start_value=0.3,
        upper_bound=0.5,
        unit="1/hr",
    ),
    # food: AUC +39 %, Cmax +76 %, tmax +1.25 hr at 20 mg (Kubitza2006)
    FitParameter(
        pid="GU__FED_f_solubility",
        lower_bound=1.0,
        start_value=4.0,
        upper_bound=20,
        unit="dimensionless",
    ),
    FitParameter(
        pid="GU__FED_f_emptying",
        lower_bound=0.05,
        start_value=0.25,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    # hepatic clearance, about 6-7 l/hr non-renal (Mueck2014)
    FitParameter(
        pid="LI__RIV2RX_k",
        lower_bound=1e-3,
        start_value=0.062,
        upper_bound=1.0,
        unit="1/min",
    ),
    # renal clearance relative to the filtration of the unbound drug, 2.4-3.1
    # l/hr (Kubitza2013, Mueck2014)
    FitParameter(
        pid="KI__RIVEX_f",
        lower_bound=1.0,
        start_value=7.9,
        upper_bound=30,
        unit="dimensionless",
    ),
    # tissue distribution, Vss about 50 l (Mueck2014)
    FitParameter(
        pid="Kp_riv",
        lower_bound=0.2,
        start_value=0.7,
        upper_bound=5.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="ftissue_riv",
        lower_bound=0.01,
        start_value=1.0,
        upper_bound=10,
        unit="l/min",
    ),
    # excretion of the metabolites: urine 30 %, feces 21 % of the dose
    # (Weinz2009); the split between bile and urine is identified, the renal
    # clearance of the metabolites (`KI__RXEX_f`) is not and is fixed
    FitParameter(
        pid="LI__RXEX_BI_k",
        lower_bound=1e-3,
        start_value=0.16,
        upper_bound=10,
        unit="1/min",
    ),
]

#: the pharmacodynamics, on the adopted pharmacokinetics
parameters_pd: list[FitParameter] = [
    # the PT is linear in the concentration: PT ratio 1.37 at 10 mg, 2.03 at
    # 80 mg (Kubitza2005a)
    FitParameter(
        pid="PT_slope",
        lower_bound=100,
        start_value=1200,
        upper_bound=10000,
        unit="1/mM",
    ),
    FitParameter(
        pid="Emax_aPTT",
        lower_bound=0.1,
        start_value=0.95,
        upper_bound=5,
        unit="dimensionless",
    ),
    FitParameter(
        pid="EC50_riv_aPTT",
        lower_bound=1e-5,
        start_value=4e-4,
        upper_bound=1e-2,
        unit="mM",
    ),
    # the inhibition of factor Xa saturates: 33 % at 10 mg, 61 % at 80 mg
    # (Kubitza2005a)
    FitParameter(
        pid="Emax_Xa",
        lower_bound=0.5,
        start_value=0.9,
        upper_bound=1.0,
        unit="dimensionless",
    ),
    FitParameter(
        pid="EC50_riv_Xa",
        lower_bound=1e-5,
        start_value=5e-4,
        upper_bound=1e-2,
        unit="mM",
    ),
]
