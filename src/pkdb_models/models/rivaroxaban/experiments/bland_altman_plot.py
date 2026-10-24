"""Bland-Altman analysis of the coagulation tests.

The data of the analysis are the sheet `Fig3b` of `bland_altman.xlsx` of the
rivaroxaban studies of pkdb_data, copied into `bland_altman_Fig3b.tsv` of the
data of the model; the statistics and the figure are written to the results of
the model.
"""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pandas.api.typing import DataFrameGroupBy
from scipy import stats
from scipy.stats import shapiro

from pkdb_models.models.rivaroxaban import DATA_PATH_RIVAROXABAN, RIVAROXABAN_PATH

#: the statistics and the figure are written here
OUTPUT_PATH = RIVAROXABAN_PATH / "results"


def shapiro_test(data: pd.Series) -> dict[str, Any]:
    """Test the differences for normality with the Shapiro-Wilk test.

    Args:
        data: the differences, without missing values.

    Returns:
        The test statistic, the p-value and whether the differences are normal.
    """
    n = len(data)

    # Perform Shapiro-Wilk test for normality
    stat, p_value = shapiro(data)
    is_normal = p_value > 0.05

    return {
        "n_points": n,
        "shapiro_stat": stat,
        "shapiro_p": p_value,
        "distribution": "normal" if is_normal else "non-normal",
        "analysis_type": "parametric" if is_normal else "non-parametric",
        "is_normal": is_normal,
    }


def save_stat(stats_df: pd.DataFrame, output_path: Path, title: str) -> None:
    """Save the statistics as `stat_<title>.csv`.

    Args:
        stats_df: the statistics.
        output_path: directory of the file, created if it does not exist.
        title: name of the analysis.
    """
    output_path.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(output_path / f"stat_{title}.csv", index=True)
    print(f"\nStatistics saved to: {output_path / f'stat_{title}.csv'}")


def calculate_statistics(diff_data: pd.Series, force: bool) -> dict[str, Any]:
    """Calculate the Bland-Altman statistics.

    Args:
        diff_data: the differences, without missing values.
        force: use the parametric statistics even if the differences are not
            normal.

    Returns:
        The statistics by name.
    """
    # Calculate basic statistics
    stats_dict = shapiro_test(diff_data)
    n = stats_dict["n_points"]

    if stats_dict["is_normal"] or force:
        # Parametric statistics
        center = np.mean(diff_data)
        sd = np.std(diff_data, ddof=1)
        stats_dict.update(
            {
                "central_tendency_name": "mean",
                "central_tendency": center,
                "standard_deviation": sd,
                "bounds_method": "mean ± 1.96 SD",
                "lower_bound": center - 1.96 * sd,
                "upper_bound": center + 1.96 * sd,
            }
        )

        # Calculate confidence intervals
        se = sd / np.sqrt(n)
        se_loas = np.sqrt(3 * sd**2 / n)
        t_interval = stats.t.interval(confidence=0.95, df=n - 1)
        ci_mean = center + np.array(t_interval) * se
        ci_upper = stats_dict["upper_bound"] + np.array(t_interval) * se_loas
        ci_lower = stats_dict["lower_bound"] + np.array(t_interval) * se_loas

        stats_dict.update(
            {
                "ci_method": "t-distribution",
                "ci_mean_lower": ci_mean[0],
                "ci_mean_upper": ci_mean[1],
                "ci_limits_lower": ci_lower,
                "ci_limits_upper": ci_upper,
            }
        )

    else:
        # Non-parametric statistics
        center = np.median(diff_data)
        stats_dict.update(
            {
                "central_tendency_name": "median",
                "central_tendency": center,
                "bounds_method": "2.5th and 97.5th percentiles",
                "lower_bound": np.percentile(diff_data, 2.5),
                "upper_bound": np.percentile(diff_data, 97.5),
            }
        )

        # Bootstrap confidence intervals
        #  bootstrap percentile method
        # does not work for CI percentiles
        bootstrap_samples = 1000
        bootstrap_medians = np.zeros(bootstrap_samples)
        for i in range(bootstrap_samples):
            sample = np.random.choice(diff_data, size=len(diff_data), replace=True)
            bootstrap_medians[i] = np.median(sample)

        ci_median = np.percentile(bootstrap_medians, [2.5, 97.5])
        stats_dict.update(
            {
                "ci_method": "bootstrap",
                "ci_mean_lower": ci_median[0],
                "ci_mean_upper": ci_median[1],
            }
        )

    # Additional statistics
    stats_dict.update(
        {
            "minimum": np.min(diff_data),
            "maximum": np.max(diff_data),
            "points_within_bounds": np.sum(
                (diff_data >= stats_dict["lower_bound"])
                & (diff_data <= stats_dict["upper_bound"])
            )
            / n
            * 100,
            "diff_mean": center,
            "in_bounds": np.sum(
                (diff_data >= stats_dict["lower_bound"])
                & (diff_data <= stats_dict["upper_bound"])
            )
            / n
            * 100,
        }
    )

    return stats_dict


def bland_altman_plot(
    df: pd.DataFrame,
    df_diff: pd.Series,
    variables: DataFrameGroupBy,
    stats_dict: dict[str, Any],
    save_fig: bool,
    output_path: Path,
    colors: list[str],
    force: bool,
    means: str,
    title: str,
) -> None:
    """Plot the differences over the means with the limits of agreement.

    Args:
        df: the data.
        df_diff: the differences, without missing values.
        variables: the data grouped by the measurement.
        stats_dict: the statistics of the differences.
        save_fig: save the figure as `figures/plot_<title>.png`.
        output_path: directory of the figure.
        colors: colors of the groups.
        force: the parametric statistics were used even if the differences are
            not normal.
        means: column of the means.
        title: name of the analysis.
    """
    # ---- Create Bland-Altman Plot ----
    # Create a Figure, which doesn't have to be square.
    fig = plt.figure(figsize=(13, 11), layout="constrained")
    # Create the main Axes.
    ax = fig.add_subplot()

    for kv, (measurement_type, group) in enumerate(variables):
        ax.scatter(
            group.x,
            group.y,
            label=measurement_type,
            color=colors[kv],
            s=80,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
            marker="o",
            zorder=3,
        )

    # Add mean / median horizontal line
    ax.axhline(
        y=stats_dict["central_tendency"],
        color="grey",
        linestyle="-",
        linewidth=2,
        zorder=2,
    )

    x_pos = df[means].max() + 0.16

    ax.text(
        x_pos,
        stats_dict["central_tendency"],
        stats_dict["central_tendency_name"].capitalize(),
        va="bottom",
        ha="right",
        fontsize=19,
        color="black",
        weight="bold",
    )

    ax.text(
        x_pos,
        stats_dict["central_tendency"],
        f"{stats_dict['central_tendency']:.1e}%",
        va="top",
        ha="right",
        fontsize=14,
        color="black",
    )

    # Add confidence intervals (limits of agreement)
    if stats_dict["is_normal"] or force:
        loa_labels = ["-1.96SD", "+1.96SD"]
    else:
        loa_labels = ["2.5th percentile", "97.5th percentile"]
    for bound, label in zip(["lower_bound", "upper_bound"], loa_labels, strict=True):
        ax.axhline(
            stats_dict[bound], linestyle="--", color="grey", linewidth=1, zorder=2
        )
        ax.text(
            x_pos,
            stats_dict[bound],
            label,
            va="bottom",
            ha="right",
            fontsize=17,
            color="black",
        )
        ax.text(
            x_pos,
            stats_dict[bound],
            f"{stats_dict[bound]:.1e}%",
            va="top",
            ha="right",
            fontsize=14,
            color="black",
        )

    # Add CIs for Mean / Median and Limits of Agreement
    # CI for the Bias (Mean / Median)
    ax.axhspan(
        stats_dict["ci_mean_lower"],
        stats_dict["ci_mean_upper"],
        color="grey",
        alpha=0.2,
    )

    if stats_dict["is_normal"] or force:
        # CI for Lower Limit of Agreement
        ax.axhspan(
            stats_dict["ci_limits_lower"][0],
            stats_dict["ci_limits_lower"][1],
            color="grey",
            alpha=0.1,
        )

        # CI for Upper Limit of Agreement
        ax.axhspan(
            stats_dict["ci_limits_upper"][0],
            stats_dict["ci_limits_upper"][1],
            color="grey",
            alpha=0.1,
        )

    # ---- Add Histogram to the Right ----
    # Use a divider to create a side-axis that scales with the main plot
    divider = make_axes_locatable(ax)
    ax_histy = divider.append_axes("right", size="15%", pad=0.2, sharey=ax)

    ax_histy.tick_params(axis="y", labelleft=False)

    # Calculate bins and weights for percentage
    y_limit = max(abs(stats_dict["lower_bound"]), abs(stats_dict["upper_bound"])) * 2.0
    # Sturges' rule
    n_bins = int(np.ceil(np.log2(len(df_diff)) + 1))
    bins = np.linspace(-y_limit, y_limit, n_bins)
    weights = np.ones_like(df_diff) * 100 / len(df_diff)  # Convert to percentages

    # Create horizontal histogram
    ax_histy.hist(
        df_diff,
        orientation="horizontal",
        weights=weights,
        color="grey",
        edgecolor="black",
        alpha=0.7,
        zorder=20,
        bins=bins,
    )

    # Add grid
    ax_histy.grid(True, axis="x", color="#EEEEEE", linewidth=0.8, zorder=1)

    # Remove y-axis ticks
    ax_histy.tick_params(axis="y", labelleft=False, left=False)

    ax_histy.set_axis_off()

    # ---- Style ----
    ax.legend(fontsize=18)
    ax.set_xlim(-0.05, 1.2)
    ax.set_ylim(-125, 125)

    # Remove the top and right spines to delete the "box"
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_title("Bland-Altman Plot", fontsize=24, pad=22, weight="bold")
    ax.set_xlabel("Mean [-]", fontsize=22, labelpad=10)  # Average of two measures
    ax.set_ylabel(
        "Difference [% mean]", fontsize=22, labelpad=10
    )  # Difference between two measures

    ax.grid(visible=True, color="#EEEEEE", linewidth=0.8, zorder=1)

    # Format labels to scientific notation
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.xaxis.set_major_formatter(formatter)

    ax.tick_params(axis="both", which="major", labelsize=14, size=5, width=2)

    ax_histy.tick_params(axis="both", which="major", labelsize=14, size=5, width=2)

    plt.tight_layout()

    if save_fig:
        fig_path = output_path / "figures"
        fig_path.mkdir(parents=True, exist_ok=True)

        fig.savefig(fig_path / f"plot_{title}.png", dpi=300)
        fig.savefig(fig_path / f"plot_{title}.svg")
    plt.close(fig)


def main() -> None:
    """Bland-Altman analysis of the data of the model."""
    data_path = DATA_PATH_RIVAROXABAN / "bland_altman_Fig3b.tsv"

    # define column with needed data
    groups = "x_measurement_type"
    differences = "y"
    means = "x"
    force_normality_calculation = True
    title = "bland_altman_mean"

    df = pd.read_csv(data_path, sep="\t")
    stats_dict = calculate_statistics(
        df[differences].dropna(),
        force=force_normality_calculation,
    )

    stats_df = pd.DataFrame.from_dict(stats_dict, orient="index", columns=["Value"])
    save_stat(stats_df, OUTPUT_PATH, title)
    # Print table to terminal
    print("\n" + "=" * 30)
    print("BLAND-ALTMAN STATISTICS")
    print("-" * 30)
    print(stats_df)
    print("=" * 30 + "\n")

    bland_altman_plot(
        df=df,
        df_diff=df[differences].dropna(),
        variables=df.groupby(groups),
        stats_dict=stats_dict,
        save_fig=True,
        output_path=OUTPUT_PATH,
        colors=["tab:purple", "tab:orange", "tab:green", "tab:red", "tab:blue"],
        force=force_normality_calculation,
        means=means,
        title=title,
    )


if __name__ == "__main__":
    main()
