# make_visuals_targets.py
# -----------------------
# Reads:
#   - dedicated_timeslot_grid_results_final.csv
#   - flexible_timeslot_grid_results_final.csv
#
# Targets (EDIT THESE):
#   - mean_avg_overtime_per_day_min <= OT_MAX
#   - mean_avg_extra_wait_min <= EX_MAX
# Objective:
#   - minimize mean_avg_wait_days (tie-break: lower overtime/day, then lower extra wait)
#
# Outputs:
#   figures_targets/
#     - dedicated_feasible_scatter.png
#     - flexible_feasible_scatter.png
#     - feasible_slack_dedicated.png
#     - feasible_slack_flexible.png
#     - overlay_ot_vs_wait.png
#     - best_vs_best_kpis.png
#   - best_under_targets.csv
#   - dedicated_feasible_sorted.csv
#   - flexible_feasible_sorted.csv

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =============================
# TARGETS (your updated values)
# =============================
OT_MAX = 0.5   # minutes/day
EX_MAX = 2.0   # minutes

# Input files
DED_PATH = "dedicated_timeslot_grid_results_final.csv"
FLEX_PATH = "flexible_timeslot_grid_results_final.csv"

REQUIRED_COLS = [
    "timeslot_type1",
    "timeslot_type2",
    "mean_avg_overtime_per_day_min",
    "mean_avg_wait_days",
    "mean_avg_extra_wait_min",
]


def load_results(path: str, policy_name: str) -> pd.DataFrame:
    df = pd.read_csv(path)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{policy_name}: missing columns {missing} in {path}")

    df = df.copy()
    df["policy"] = policy_name
    df["feasible"] = (
        (df["mean_avg_overtime_per_day_min"] <= OT_MAX) &
        (df["mean_avg_extra_wait_min"] <= EX_MAX)
    )
    return df


def sort_feasible(df: pd.DataFrame) -> pd.DataFrame:
    feas = df[df["feasible"]].copy()
    feas = feas.sort_values(
        by=[
            "mean_avg_wait_days",                 # objective
            "mean_avg_overtime_per_day_min",      # tie-break 1
            "mean_avg_extra_wait_min",            # tie-break 2
        ],
        ascending=True
    ).reset_index(drop=True)
    return feas


def best_row(feas_sorted: pd.DataFrame):
    if len(feas_sorted) == 0:
        return None
    return feas_sorted.iloc[0]


def _prep_style():
    # Light, clean grid without forcing weird fonts/colors
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.alpha"] = 0.25


def _annotate_best(best: pd.Series, x_col: str, y_col: str):
    bx = float(best[x_col])
    by = float(best[y_col])
    ts1 = int(best["timeslot_type1"])
    ts2 = int(best["timeslot_type2"])
    plt.scatter([bx], [by], marker="*", s=260, edgecolors="black", linewidths=0.6, zorder=5)
    plt.text(bx, by, f"  best=({ts1},{ts2})", fontsize=10, va="center")


def scatter_ot_vs_wait(df: pd.DataFrame, best: pd.Series, title: str, outpath: str) -> None:
    x = df["mean_avg_overtime_per_day_min"].to_numpy()
    y = df["mean_avg_wait_days"].to_numpy()
    ex = df["mean_avg_extra_wait_min"].to_numpy()
    feasible = df["feasible"].to_numpy()

    plt.figure()

    # Infeasible: faint gray x
    if np.any(~feasible):
        plt.scatter(
            x[~feasible], y[~feasible],
            marker="x", s=40, alpha=0.35,
            label="Infeasible"
        )

    # Feasible: colored by EX with colorbar
    if np.any(feasible):
        sc = plt.scatter(
            x[feasible], y[feasible],
            marker="o", s=70, alpha=0.85,
            c=ex[feasible], cmap="viridis",
            edgecolors="white", linewidths=0.4,
            label="Feasible (colored by EX)"
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("mean_avg_extra_wait_min")

    # Constraint line for OT
    plt.axvline(OT_MAX, linestyle="--", linewidth=1.5, alpha=0.9, label=f"OT target = {OT_MAX}")

    plt.title(title)
    plt.xlabel("mean_avg_overtime_per_day_min")
    plt.ylabel("mean_avg_wait_days")
    plt.legend(loc="best")
    plt.tight_layout()

    if best is not None:
        _annotate_best(best, "mean_avg_overtime_per_day_min", "mean_avg_wait_days")
        plt.tight_layout()

    plt.savefig(outpath, dpi=200)
    plt.close()


def slack_plot_ex_vs_wait(df: pd.DataFrame, best: pd.Series, title: str, outpath: str) -> None:
    # EX vs Wait, show feasible in color by OT (nice: shows the tradeoff)
    ex = df["mean_avg_extra_wait_min"].to_numpy()
    w = df["mean_avg_wait_days"].to_numpy()
    ot = df["mean_avg_overtime_per_day_min"].to_numpy()
    feasible = df["feasible"].to_numpy()

    plt.figure()

    # Infeasible: faint gray x
    if np.any(~feasible):
        plt.scatter(
            ex[~feasible], w[~feasible],
            marker="x", s=40, alpha=0.35,
            label="Infeasible"
        )

    # Feasible: colored by OT
    if np.any(feasible):
        sc = plt.scatter(
            ex[feasible], w[feasible],
            marker="o", s=70, alpha=0.85,
            c=ot[feasible], cmap="plasma",
            edgecolors="white", linewidths=0.4,
            label="Feasible (colored by OT)"
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("mean_avg_overtime_per_day_min")

    plt.axvline(EX_MAX, linestyle="--", linewidth=1.5, alpha=0.9, label=f"EX target = {EX_MAX}")

    plt.title(title)
    plt.xlabel("mean_avg_extra_wait_min")
    plt.ylabel("mean_avg_wait_days")
    plt.legend(loc="best")
    plt.tight_layout()

    if best is not None:
        _annotate_best(best, "mean_avg_extra_wait_min", "mean_avg_wait_days")
        plt.tight_layout()

    plt.savefig(outpath, dpi=200)
    plt.close()


def overlay_plot(ded: pd.DataFrame, flex: pd.DataFrame, outpath: str) -> None:
    # Overlay feasible points from both policies, color by EX
    plt.figure()

    for df, label, marker in [
        (ded, "Dedicated", "o"),
        (flex, "Flexible", "s"),
    ]:
        feasible = df["feasible"].to_numpy()
        if not np.any(feasible):
            continue

        sc = plt.scatter(
            df.loc[feasible, "mean_avg_overtime_per_day_min"],
            df.loc[feasible, "mean_avg_wait_days"],
            c=df.loc[feasible, "mean_avg_extra_wait_min"],
            cmap="viridis",
            marker=marker,
            s=70,
            alpha=0.85,
            edgecolors="white",
            linewidths=0.4,
            label=f"{label} feasible"
        )

    # One colorbar for EX
    cbar = plt.colorbar(sc)
    cbar.set_label("mean_avg_extra_wait_min")

    plt.axvline(OT_MAX, linestyle="--", linewidth=1.5, alpha=0.9, label=f"OT target = {OT_MAX}")

    plt.title(f"Feasible region overlay (OT ≤ {OT_MAX}, EX ≤ {EX_MAX})")
    plt.xlabel("mean_avg_overtime_per_day_min")
    plt.ylabel("mean_avg_wait_days")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def plot_best_vs_best(ded_best: pd.Series, flex_best: pd.Series, outpath: str) -> None:
    labels = ["Avg wait (days)", "OT/day (min)", "EX (min)"]

    ded_vals = [
        float(ded_best["mean_avg_wait_days"]),
        float(ded_best["mean_avg_overtime_per_day_min"]),
        float(ded_best["mean_avg_extra_wait_min"]),
    ]
    flex_vals = [
        float(flex_best["mean_avg_wait_days"]),
        float(flex_best["mean_avg_overtime_per_day_min"]),
        float(flex_best["mean_avg_extra_wait_min"]),
    ]

    x = np.arange(len(labels))
    width = 0.38

    plt.figure(figsize=(9, 5))
    plt.bar(x - width / 2, ded_vals, width, label="Dedicated", alpha=0.85)
    plt.bar(x + width / 2, flex_vals, width, label="Flexible", alpha=0.85)

    # Annotate values
    for i, v in enumerate(ded_vals):
        plt.text(i - width / 2, v, f"  {v:.3g}", va="bottom", fontsize=10)
    for i, v in enumerate(flex_vals):
        plt.text(i + width / 2, v, f"  {v:.3g}", va="bottom", fontsize=10)

    plt.xticks(x, labels)
    plt.title(f"Best under targets (OT ≤ {OT_MAX}, EX ≤ {EX_MAX}), objective: minimize waiting time")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


if __name__ == "__main__":
    _prep_style()

    outdir = "figures_targets"
    os.makedirs(outdir, exist_ok=True)

    # ---- load ----
    ded = load_results(DED_PATH, "dedicated")
    flex = load_results(FLEX_PATH, "flexible")

    # ---- feasible sets + best ----
    ded_feas = sort_feasible(ded)
    flex_feas = sort_feasible(flex)

    ded_best = best_row(ded_feas)
    flex_best = best_row(flex_feas)

    # Save feasible tables (sorted by objective/tie-breaks)
    ded_feas.to_csv("dedicated_feasible_sorted.csv", index=False)
    flex_feas.to_csv("flexible_feasible_sorted.csv", index=False)

    # Summary CSV
    summary_rows = []
    for name, df_all, df_feas, best in [
        ("dedicated", ded, ded_feas, ded_best),
        ("flexible", flex, flex_feas, flex_best),
    ]:
        summary_rows.append({
            "policy": name,
            "num_total_combos": int(len(df_all)),
            "num_feasible_combos": int(len(df_feas)),
            "best_timeslot_type1": int(best["timeslot_type1"]) if best is not None else None,
            "best_timeslot_type2": int(best["timeslot_type2"]) if best is not None else None,
            "best_mean_avg_wait_days": float(best["mean_avg_wait_days"]) if best is not None else None,
            "best_mean_avg_overtime_per_day_min": float(best["mean_avg_overtime_per_day_min"]) if best is not None else None,
            "best_mean_avg_extra_wait_min": float(best["mean_avg_extra_wait_min"]) if best is not None else None,
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv("best_under_targets.csv", index=False)

    print(f"\nTargets: OT/day <= {OT_MAX} min | EX <= {EX_MAX} min | Objective: minimize mean_avg_wait_days")
    print(summary)

    # ---- Figures: OT vs Wait (feasible highlighted, colored by EX) ----
    scatter_ot_vs_wait(
        ded,
        ded_best,
        title=f"Dedicated: OT/day vs Waiting (feasible: OT≤{OT_MAX}, EX≤{EX_MAX}; color=EX)",
        outpath=os.path.join(outdir, "dedicated_feasible_scatter.png"),
    )

    scatter_ot_vs_wait(
        flex,
        flex_best,
        title=f"Flexible: OT/day vs Waiting (feasible: OT≤{OT_MAX}, EX≤{EX_MAX}; color=EX)",
        outpath=os.path.join(outdir, "flexible_feasible_scatter.png"),
    )

    # ---- Figures: Slack plot EX vs Wait (feasible colored by OT) ----
    slack_plot_ex_vs_wait(
        ded,
        ded_best,
        title=f"Dedicated: EX vs Waiting (feasible: OT≤{OT_MAX}, EX≤{EX_MAX}; color=OT)",
        outpath=os.path.join(outdir, "feasible_slack_dedicated.png"),
    )

    slack_plot_ex_vs_wait(
        flex,
        flex_best,
        title=f"Flexible: EX vs Waiting (feasible: OT≤{OT_MAX}, EX≤{EX_MAX}; color=OT)",
        outpath=os.path.join(outdir, "feasible_slack_flexible.png"),
    )

    # ---- Figure: overlay feasible points ----
    overlay_plot(
        ded,
        flex,
        outpath=os.path.join(outdir, "overlay_ot_vs_wait.png"),
    )

    # ---- Figure: best vs best ----
    if ded_best is not None and flex_best is not None:
        plot_best_vs_best(
            ded_best,
            flex_best,
            outpath=os.path.join(outdir, "best_vs_best_kpis.png"),
        )
    else:
        print("\nNote: best_vs_best_kpis.png not created because one policy has 0 feasible solutions under these targets.")

    print("\nSaved figures to:", outdir)
    print("Saved tables: best_under_targets.csv, dedicated_feasible_sorted.csv, flexible_feasible_sorted.csv")
