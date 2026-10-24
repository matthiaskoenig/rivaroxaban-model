"""Parameter fit problems for rivaroxaban.

The data of a fit is selected from the fit mappings of the studies in three
steps, see `sbmlsim.fit.helpers`: the filters on the metadata of the mappings
select the training data of the fit and exclude the rest, the outliers are
named once for all mappings and are the training data which is not usable, and
the validation data is the part of the selected data a fit is evaluated on but
not fitted to. The metadata of a mapping describes its curve, what a fit does
with the curve is decided here.

The training data are the single doses of healthy subjects without an
interacting drug, fasted or fed: the oral solutions separate the distribution
and the clearance from the dissolution, the tablets of 1.25 to 80 mg the
saturable dissolution, the fed arms the effects of food, the urine the renal
clearance and the mass balance of [14C]rivaroxaban (Weinz2009) the excretion of
the metabolites. The validation data has four parts, each a different
question: the interactions (`INTERACTION_VALIDATION`, the arms with a
coadministered drug whose interaction the model describes with factors which
are not fitted), the doses (`DOSE_VALIDATION`, the repeated doses, which the
single doses of the training predict), the populations
(`POPULATION_VALIDATION`, renal and mild hepatic impairment, which the model
describes by the renal function and the cirrhosis) and the limitations
(`LIMITATION_VALIDATION`, data of a population or a measurement the model does
not describe). The summaries of data which is fitted one by one and the second
curve of the same subjects are excluded (`AGGREGATE_MAPPINGS`).
"""

from pkpdsim.fit import metadata_filter, observable_filter
from pkpdsim.metadata import Coadministration, Dosing, Health, Route
from sbmlsim.console import console
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.fit import FitMapping, FitMappingCollection, MappingKind
from sbmlsim.fit.helpers import FitMappings, MappingFilter

from pkdb_models.models.rivaroxaban import DATA_PATHS, RIVAROXABAN_PATH
from pkdb_models.models.rivaroxaban.experiments.studies import (
    Alalawneh2024,
    Ding2020,
    Frost2014,
    GenisNajera2022,
    Graff2007,
    Greenblatt2018,
    Jiang2010,
    Kreutz2017,
    Kubitza2005,
    Kubitza2005a,
    Kubitza2006,
    Kubitza2013,
    Lenard2024,
    Lenard2025,
    Liu2024,
    Metze2021,
    Moore2014,
    Rohr2024,
    Santos2024,
    Tang2025,
    Weinz2009,
    Zhao2009,
)

#: the keys of fit mappings by study: the keys of the mappings are SIds of their
#: experiment, which two studies can share, so a named mapping is named with
#: its study
type NamedMappings = dict[str, set[str]]

EXPERIMENT_CLASSES: list[type[SimulationExperiment]] = [
    Alalawneh2024,
    Ding2020,
    Frost2014,  # the rivaroxaban groups of the apixaban study
    GenisNajera2022,
    Graff2007,
    Greenblatt2018,
    Jiang2010,
    Kreutz2017,  # the rivaroxaban groups of the apixaban study
    Kubitza2005,
    Kubitza2005a,
    Kubitza2006,
    Kubitza2013,
    Lenard2024,  # the rivaroxaban groups of the edoxaban study
    Lenard2025,  # the rivaroxaban groups of the edoxaban study
    Liu2024,
    Metze2021,  # the rivaroxaban group of the edoxaban study
    Moore2014,
    Rohr2024,  # the rivaroxaban groups of the edoxaban study
    Santos2024,
    Tang2025,
    Weinz2009,
    Zhao2009,
]

#: observables of the pharmacokinetics
PK_OBSERVABLES: set[str] = {
    "Cve_riv",
    "Cve_rx",
    "Cve_riv_rx",
    "Aurine_riv",
    "Aurine_rx",
    "Aurine_riv_rx",
    "Afeces_riv",
    "Afeces_rx",
    "Afeces_riv_rx",
    "Afeces_urine_riv_rx",
    "KI__RIVEX",
}

#: observables of the pharmacodynamics
PD_OBSERVABLES: set[str] = {
    "PT",
    "PT_change",
    "PT_ratio",
    "aPTT",
    "aPTT_change",
    "aPTT_ratio",
    "Xa_inhibition",
}

#: mappings which are not fitted, the data is not usable
OUTLIER_MAPPINGS: NamedMappings = {
    # the inhibition of factor Xa of the placebo groups: the model is exactly
    # 0 without a dose, whatever its parameters, the data scatter by -6 to 15 %
    # around 0 with small SDs, which made 90 % of the cost of the PD fit
    # (PD_20261018_082110__d163a) without information on a parameter
    "Jiang2010": {"fm_Fig2_factor_Xa_inhibition_P"},
    "Kubitza2005": {"fm_Fig1A_factor_Xa_inhibition_P"},
    "Kubitza2005a": {"fm_Fig3A_factor_Xa_inhibition_P"},
    "Zhao2009": {"fm_factor_Xa_inhibition_SDP", "fm_factor_Xa_inhibition_MDP"},
}

#: the validation of the interactions: the arms with a coadministered drug,
#: simulated with literature based factors of the hepatic CYP3A4
#: (`LI__f_cyp3a4`) and the renal P-glycoprotein (`KI__f_pgp`) which are not
#: fitted. The model has no intestinal CYP3A4 and P-glycoprotein, so the rise
#: of Cmax of the inhibitors is not described
INTERACTION_VALIDATION: NamedMappings = {
    # verapamil ER 360 mg at steady state, CYP3A4 0.3 (AUC 1.39-1.43)
    "Greenblatt2018": {
        "fm_po_NRF_VER_rivaroxaban",
        "fm_po_NRF_VER_rivaroxaban_urine",
        "fm_po_MRI_VER_rivaroxaban",
        "fm_po_MRI_VER_rivaroxaban_urine",
    },
    # clarithromycin 500 mg bid, CYP3A4 0.1, P-gp 0.7 (AUC 1.43)
    "Lenard2024": {"fm_riv25_clar"},
    # carbamazepine, induction of CYP3A4 and P-gp (AUC 0.56); gabapentin and
    # pregabalin, negative controls without an interaction
    "Lenard2025": {"fm_riv25_car", "fm_riv25_gab", "fm_riv25_pre"},
    # erythromycin 500 mg tid, CYP3A4 0.3, P-gp 0.7 (AUC 1.38 with normal
    # renal function, 1.60 and 1.53 in renal impairment)
    "Moore2014": {
        "fm_rivaroxaban_NRF_ERY",
        "fm_rivaroxaban_mildRI_ERY",
        "fm_rivaroxaban_moderateRI_ERY",
        "fm_rivaroxaban_urine_NRF_ERY",
        "fm_rivaroxaban_urine_mildRI_ERY",
        "fm_rivaroxaban_urine_moderateRI_ERY",
    },
    # ritonavir 100 mg bid, CYP3A4 0.05, P-gp 0.65 (AUC 1.87, CLrenal 0.70)
    "Rohr2024": {
        "fm_rivaroxaban_API25__RIV25__EDO50__RITx5",
        "fm_rivaroxaban_urine_API25__RIV25__EDO50__RITx5",
    },
}

#: the validation of the doses: the repeated doses, which the single doses of
#: the training predict
DOSE_VALIDATION: NamedMappings = {
    # 10 mg once daily, the last dose
    "Frost2014": {"fm_rivaroxaban_mRIV10"},
    # 20 mg once daily for 7 days, fed
    "Kreutz2017": {
        "fm_PO_riv_rivaroxaban",
        "fm_PO_riv_prothrombin_time_ratio",
        "fm_PO_riv_aPTT_ratio",
    },
    # 5-30 mg twice daily and 5 mg once and three times daily for 7 days, fed
    "Kubitza2005": {
        *(f"fm_Fig3_rivaroxaban_BID{dose}" for dose in [5, 10, 20, 30]),
        *(
            f"fm_Fig1{panel}_{sid}_{group}"
            for panel, sid in [
                ("A", "factor_Xa_inhibition"),
                ("B", "prothrombin_time__change_relative_"),
                ("C", "aPTT__change_relative_"),
            ]
            for group in ["P", "BID5", "BID10", "BID20", "BID30"]
            # the inhibition of factor Xa of the placebo is an outlier
            if (panel, group) != ("A", "P")
        ),
        *(
            f"fm_Fig2{panel}_{sid}_{group}"
            for panel, sid in [
                ("A", "factor_Xa_inhibition"),
                ("B", "prothrombin_time__change_relative_"),
                ("C", "aPTT__change_relative_"),
            ]
            for group in ["OD5", "TID5"]
        ),
    },
    # 5-30 mg for 6 days (10 doses), fed
    "Zhao2009": {
        *(
            f"fm_{sid}_MD{dose}"
            for sid in ["rivaroxaban", "rivaroxaban_urine"]
            for dose in [5, 10, 20, 30]
        ),
        *(
            f"fm_{sid}_MD{group}"
            for sid in ["factor_Xa_inhibition", "prothrombin_time__change_relative_"]
            for group in ["P", "5", "10", "20", "30"]
            # the inhibition of factor Xa of the placebo is an outlier
            if (sid, group) != ("factor_Xa_inhibition", "P")
        ),
    },
}

#: the validation of the populations: renal impairment (the renal function
#: relative to a GFR of 101 ml/min) and mild hepatic impairment (Child-Pugh A,
#: the shared cirrhosis model)
POPULATION_VALIDATION: NamedMappings = {
    # moderate renal impairment, CrCl 71 ml/min
    "Greenblatt2018": {"fm_po_MRI_rivaroxaban", "fm_po_MRI_rivaroxaban_urine"},
    "Kubitza2013": {
        "fm_Fig1_rivaroxaban_mildHI",
        "fm_Fig2A_factor_Xa_inhibition_mildHI",
        "fm_Fig2B_prothrombin_time__change_relative__mildHI",
    },
    # elderly patients with a CrCl of 67 ml/min, anti-Xa calibrated levels
    "Metze2021": {"fm_rivaroxaban_RIV20"},
    # mild and moderate renal impairment, CrCl 67.1 and 42.5 ml/min
    "Moore2014": {
        "fm_rivaroxaban_mildRI",
        "fm_rivaroxaban_moderateRI",
        "fm_rivaroxaban_urine_mildRI",
        "fm_rivaroxaban_urine_moderateRI",
    },
}

#: the validation of the limitations of the model: valid data of a population
#: or a measurement the model does not describe, which would distort the fitted
#: parameters of the processes it does describe
LIMITATION_VALIDATION: NamedMappings = {
    "Alalawneh2024": {
        # the obese: the exposure is independent of the body weight, the
        # model scales the organs and the hepatic clearance with it
        "fm_rivaroxaban_obese",
        "fm_rivaroxaban_urine_obese",
        "fm_rivaroxaban_excretion_urine_obese",
        # the absolute PT and aPTT of the laboratory without SD, the noise of
        # the aPTT (5 s) is larger than the effect of the drug after 18 hr
        "fm_prothrombin_time_non_obese",
        "fm_aPTT_non_obese",
        "fm_prothrombin_time_obese",
        "fm_aPTT_obese",
    },
    # moderate hepatic impairment (Child-Pugh B): AUC x1.80 vs 2.27, the
    # urinary fraction rises in the model while the renal clearance falls
    # (2.4 to 0.7 l/hr), the PD baseline reflects the reduced synthesis of the
    # coagulation factors
    "Kubitza2013": {
        "fm_Fig1_rivaroxaban_moderateHI",
        "fm_Fig2A_factor_Xa_inhibition_moderateHI",
        "fm_Fig2B_prothrombin_time__change_relative__moderateHI",
    },
}

#: the summaries of data which the fits use one by one and the second curve of
#: the same subjects, excluded from every fit: they would count the subjects
#: twice
AGGREGATE_MAPPINGS: NamedMappings = {
    # the morbidly obese are 6 of the 18 obese; the excretion rates are the
    # interval amounts of the cumulative urine
    "Alalawneh2024": {
        "fm_rivaroxaban_morbidly_obese",
        "fm_rivaroxaban_excretion_urine_non_obese",
    },
    # the test tablets of the replicate crossover, the same subjects as the
    # reference tablets
    "Ding2020": {"fm_po_Fasting_T10", "fm_po_Fed_T10"},
    "GenisNajera2022": {"fm_po_T20"},
    "Santos2024": {"fm_po_Fasting_T20", "fm_po_Fed_T20"},
    "Tang2025": {
        "fm_po_Fasting_T1",
        "fm_po_Fasting_R2",
        "fm_po_Fasting_T2",
        "fm_po_Fed_T1",
        "fm_po_Fed_R2",
        "fm_po_Fed_T2",
    },
    # the 20 patients of the mean
    "Metze2021": {f"fm_rivaroxaban_RIV20_IS{k:03d}" for k in range(41, 61)},
    # the sum of the urine and the feces
    "Weinz2009": {"fm_riv_rx_feces_urine"},
}

#: the validation data
VALIDATION_MAPPINGS: NamedMappings = {
    study: set().union(
        *(
            mappings.get(study, set())
            for mappings in [
                INTERACTION_VALIDATION,
                DOSE_VALIDATION,
                POPULATION_VALIDATION,
                LIMITATION_VALIDATION,
            ]
        )
    )
    for study in INTERACTION_VALIDATION
    | DOSE_VALIDATION
    | POPULATION_VALIDATION
    | LIMITATION_VALIDATION
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


#: the training data: single doses of oral rivaroxaban in healthy subjects
#: without a coadministered drug, fasted or fed
filter_training = metadata_filter(
    route={Route.PO, Route.IV},
    coadministration={Coadministration.NONE},
    health={Health.HEALTHY},
    excluded={"dosing": {Dosing.MULTIPLE}},
)


def filter_selection(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """The training data and the named validation data, nothing else."""
    if named(VALIDATION_MAPPINGS, fit_mapping_key, fit_mapping):
        return True
    return filter_training(fit_mapping_key, fit_mapping)


def filter_aggregates(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """No summaries of data which is fitted one by one, see `AGGREGATE_MAPPINGS`."""
    return not named(AGGREGATE_MAPPINGS, fit_mapping_key, fit_mapping)


#: Only the pharmacokinetic data.
filter_pk = observable_filter(PK_OBSERVABLES)


#: Only the pharmacodynamic data.
filter_pd = observable_filter(PD_OBSERVABLES)


def mapping_collections(
    filters: MappingFilter | list[MappingFilter],
    validation: NamedMappings = VALIDATION_MAPPINGS,
) -> dict[str, list[FitMappingCollection]]:
    """Select the data of a fit from the fit mappings of the studies.

    The filters select the training data, the mappings which fail them are
    excluded, and so are the summaries of data which is fitted one by one
    (`AGGREGATE_MAPPINGS`) and every mapping which is neither training data
    (`filter_training`) nor named validation data (`filter_selection`). The
    outliers are `OUTLIER_MAPPINGS`, the validation data is
    `VALIDATION_MAPPINGS`. Neither kind enters the cost; both are simulated and
    evaluated with the training data, so a report says how the fitted
    parameters describe data they were not fitted on.

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
        base_path=RIVAROXABAN_PATH,
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
            filter_selection,
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
                ("populations", POPULATION_VALIDATION),
                ("limitations", LIMITATION_VALIDATION),
            ]
        )
    )
    return selection.collections


def f_collections_pk() -> dict[str, list[FitMappingCollection]]:
    """Pharmacokinetic data."""
    return mapping_collections([filter_pk])


def f_collections_pd() -> dict[str, list[FitMappingCollection]]:
    """Pharmacodynamic data."""
    return mapping_collections([filter_pd])


if __name__ == "__main__":
    for f in [f_collections_pk, f_collections_pd]:
        console.rule(title=f.__name__, align="left", style="white")
        f()
