"""Pharmacokinetic and pharmacodynamic parameters of the rivaroxaban scans.

The pharmacokinetic parameters are the non-compartmental analysis of
`pkpdutils`, named like the variables of its `NCAResult` (e.g. `auc_last`,
`auc_inf_obs`, `cmax`, `thalf`, `lambda_z`, `cl_f`, `vz_f`), with the unit of
every parameter in `<name>_unit`, for rivaroxaban, its metabolites and their
sum in plasma. The dose of all three is the oral dose of rivaroxaban (tablet
and solution). The concentrations below `LLOQ_FRACTION` of the maximum are
below the limit of quantification, so that the noise of the integrator is not
the terminal phase. The pharmacodynamic parameters are the minimum and the
maximum of the effects.

A scan has one dimension, the scanned parameter; a run over the draws of the
fitted parameters (`run_simulation pkdb_models.models.rivaroxaban
--uncertainty`) adds the dimension of the draws. The parameters are calculated
for every point of every dimension, a row per point with its position along
each dimension (the columns `point` and, with draws, `draw`), see `summarize`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import xarray as xr
from pkpdsim.pk import nca_frame, scan_nca
from pkpdsim.scan import end_values, scan_dimensions, start_values
from pkpdutils import Route, ureg
from sbmlsim import Q
from sbmlsim.result import ScanResult

if TYPE_CHECKING:
    from pkdb_models.models.rivaroxaban.experiments.base_experiment import (
        RivaroxabanSimulationExperiment,
    )

#: the oral doses of the scans, of a tablet and of a solution
DOSES: tuple[str, ...] = ("PODOSE_riv", "POSOLUTION_riv")
#: the substances by id of the concentration in plasma (`[Cve_<id>]`)
SUBSTANCES: dict[str, str] = {
    "riv": "rivaroxaban",
    "rx": "rivaroxaban metabolites",
    "riv_rx": "rivaroxaban and metabolites",
}
#: the effects whose minimum and maximum are the pharmacodynamic parameters
EFFECTS: list[str] = [
    "PT",
    "aPTT",
    "Xa_inhibition",
]
#: the limit of quantification relative to the maximal concentration
LLOQ_FRACTION = 1e-4


def _rows(
    df: pd.DataFrame, res: ScanResult, across: str | None, **columns: str
) -> pd.DataFrame:
    """Get the rows of the parameters of a scan, in the order of the points.

    The scan has one dimension besides the draws `across`; its position is the
    column `point`, the position of the draw the column `draw`, followed by
    the `columns` (e.g. the substance); the rows are ordered by the point and
    then the draw.

    Raises:
        ValueError: if the scan has not one dimension besides the draws.
    """
    dims = scan_dimensions(res, across)
    if len(dims) != 1:
        raise ValueError(
            f"A scan of one dimension is required, not of {dims} (draws: {across})."
        )
    names = {dims[0]: "point", **({across: "draw"} if across is not None else {})}
    df = df.rename(columns=names)
    positions = list(names.values())
    df = df.sort_values(positions, kind="stable", ignore_index=True)
    for k, (name, value) in enumerate(columns.items()):
        df.insert(k, name, value)
    return df[[*positions, *(c for c in df.columns if c not in positions)]]


def pk_parameters(
    experiment: RivaroxabanSimulationExperiment,
    res: ScanResult,
    across: str | None,
    sid: str,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of a substance in every simulation of a scan.

    The non-compartmental analysis of all simulations at once. The renal and
    fecal clearances are the amounts excreted in the urine and the feces until
    the end of the simulation over `auc_last`, the total clearance by
    excretion is their sum.

    Args:
        experiment: experiment of the scan, for the molecular weight.
        res: result of the scan.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.
        sid: id of the substance, a key of `SUBSTANCES`.

    Returns:
        One row per simulation with its position (`point`, and `draw` with
        draws), the id of the substance (`substance`), the parameters and their
        units (`<name>_unit`).
    """
    units = res.units

    # the dose in the substance unit of the concentration, which `cl_f` and
    # `vz_f` require
    dose_mass = sum(
        (Q(start_values(res, key).values, units[key]) for key in DOSES if key in units),
        start=Q(0, "mg"),
    )
    result = scan_nca(
        res,
        f"[Cve_{sid}]",
        dose=dose_mass / experiment.Mr.riv,
        route=Route.ORAL,
        substance=SUBSTANCES[sid],
        lloq_fraction=LLOQ_FRACTION,
    )
    df = nca_frame(result).drop(columns="nca_flags")
    position = tuple(df[dim].to_numpy() for dim in res.dims)

    excretions = {"cl_renal": f"Aurine_{sid}", "cl_fecal": f"Afeces_{sid}"}
    amounts = {
        amount: ureg.Quantity(
            end_values(res, amount).values[position], str(units[amount])
        )
        for amount in excretions.values()
    }
    auc_last = ureg.Quantity(df["auc_last"].to_numpy(), df["auc_last_unit"].iloc[0])
    clearances = {
        cl_key: (amounts[amount] / auc_last).to("ml/min")
        for cl_key, amount in excretions.items()
    }
    clearances["cl_total"] = clearances["cl_renal"] + clearances["cl_fecal"]
    for key, value in {**amounts, **clearances}.items():
        df[key] = value.magnitude
        df[f"{key}_unit"] = str(value.units)
    return _rows(df, res, across, substance=sid)


def calculate_rivaroxaban_pk(
    experiment: RivaroxabanSimulationExperiment,
    res: ScanResult,
    across: str | None = None,
) -> pd.DataFrame:
    """Calculate the pharmacokinetic parameters of every simulation of a scan.

    Args:
        experiment: experiment of the scan.
        res: result of the scan.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        One row per substance and simulation of the scan, in the order of the
        simulations for every substance, with its position (`point`, and
        `draw` with draws).
    """
    return pd.concat(
        [pk_parameters(experiment, res, across, sid) for sid in SUBSTANCES],
        ignore_index=True,
    )


def calculate_rivaroxaban_pd(
    res: ScanResult, across: str | None = None
) -> dict[str, pd.DataFrame]:
    """Calculate the pharmacodynamic parameters of every simulation of a scan.

    Args:
        res: result of the scan.
        across: the dimension of the draws of the fitted parameters, `None`
            for a scan without draws.

    Returns:
        By effect, one row per simulation with its position (`point`, and
        `draw` with draws), the minimum and the maximum.
    """
    # the extrema of the timecourses of all simulations at once
    units = res.units
    frames: dict[str, pd.DataFrame] = {}
    for sid in EFFECTS:
        timecourses = res.to_timecourses(sid)
        values = timecourses.ds["value"]
        times = [d for d in values.dims if d not in timecourses.sample_dims]
        extrema = xr.Dataset(
            {"min": values.min(dim=times), "max": values.max(dim=times)}
        )
        df = extrema.to_dataframe().reset_index()
        df = _rows(df[[*res.dims, "min", "max"]], res, across, sid=sid)
        df["unit"] = str(units[sid])
        frames[sid] = df
    return frames


def summarize(
    df: pd.DataFrame, key: str, quantiles: tuple[float, float]
) -> pd.DataFrame:
    """Get the median and the quantiles of a parameter over the draws of every point.

    Args:
        df: the parameters with the columns `point` and `draw`.
        key: the parameter.
        quantiles: the lower and the upper quantile.

    Returns:
        One row per point (the index), the columns `lower`, `median`, `upper`.
    """
    low, high = quantiles
    grouped = df.groupby("point")[key]
    return pd.DataFrame(
        {
            "lower": grouped.quantile(low),
            "median": grouped.median(),
            "upper": grouped.quantile(high),
        }
    ).astype(np.float64)
