#!/usr/bin/env python3
"""
Reproduce the Monte Carlo robustness artifacts used in the IEEE Access submission.

This script regenerates:
  - Fig. 8 (figs/fig8_mc_gain.png): Monte Carlo gain intervals under finite-shot resampling
  - supplement/mc_gain_summary.csv
  - supplement/mc_coefficients_summary.csv

Inputs (in supplement/):
  - table_budget_full.csv
  - table_coefficients.csv

The random seed is fixed to make results deterministic.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SUPP = ROOT / "supplement"
FIGS = ROOT / "figs"

SEED = 20260223
SHOTS = 8192
N_REP = 30000  # number of Monte Carlo trials for gain intervals


def H2(e: np.ndarray) -> np.ndarray:
    e = np.clip(e, 1e-12, 1 - 1e-12)
    return -(e * np.log2(e) + (1 - e) * np.log2(1 - e))


def mi_from_e(e: np.ndarray) -> np.ndarray:
    return 1 - H2(e)


def e_from_mi(mi: float) -> float:
    # solve mi = 1 - H2(e), e in [0, 0.5] via bisection
    target = 1 - float(mi)
    if target <= 0:
        return 0.0
    if target >= 1:
        return 0.5
    lo, hi = 0.0, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if H2(np.array([mid]))[0] > target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def parse_baseline(s: str) -> tuple[int, int, float]:
    m = re.match(r"\((\d+)\s*,\s*(\d+)\)\s*,\s*([0-9.]+)", s.strip())
    if not m:
        raise ValueError(f"Cannot parse baseline field: {s}")
    return int(m.group(1)), int(m.group(2)), float(m.group(3))


def parse_control(s: str) -> tuple[int, int, float, float]:
    m = re.match(r"\((\d+)\s*,\s*(\d+)\s*,\s*([0-9.]+)\)\s*,\s*([0-9.]+)", s.strip())
    if not m:
        raise ValueError(f"Cannot parse control field: {s}")
    return int(m.group(1)), int(m.group(2)), float(m.group(3)), float(m.group(4))


def mc_gain_intervals() -> pd.DataFrame:
    df = pd.read_csv(SUPP / "table_budget_full.csv")
    rng = np.random.default_rng(SEED)

    rows = []
    for _, r in df.iterrows():
        p = float(r["p"])
        B = int(r["B"])
        _, _, I0 = parse_baseline(r["Baseline (n,L), I"])
        _, _, _, I1 = parse_control(r["Controlled (n,L,u), I"])

        e0 = e_from_mi(I0)
        e1 = e_from_mi(I1)

        k0 = rng.binomial(SHOTS, e0, size=N_REP)
        k1 = rng.binomial(SHOTS, e1, size=N_REP)

        I0_hat = mi_from_e(np.minimum(k0 / SHOTS, 0.5))
        I1_hat = mi_from_e(np.minimum(k1 / SHOTS, 0.5))

        gains = (I1_hat - I0_hat) / I0_hat * 100.0

        rows.append(
            dict(
                p=p,
                B=B,
                gain_median=float(np.median(gains)),
                gain_p05=float(np.quantile(gains, 0.05)),
                gain_p95=float(np.quantile(gains, 0.95)),
                gain_mean=float(np.mean(gains)),
                gain_std=float(np.std(gains, ddof=1)),
            )
        )

    out = pd.DataFrame(rows).sort_values(["p", "B"])
    out.to_csv(SUPP / "mc_gain_summary.csv", index=False)
    return out


def plot_gain(df: pd.DataFrame) -> None:
    ps = sorted(df["p"].unique())
    plt.figure(figsize=(6.2, 3.6))
    for p in ps:
        sub = df[df["p"] == p].sort_values("B")
        x = sub["B"].to_numpy()
        y = sub["gain_median"].to_numpy()
        yerr = np.vstack([y - sub["gain_p05"].to_numpy(), sub["gain_p95"].to_numpy() - y])
        plt.errorbar(x, y, yerr=yerr, marker="o", linestyle="-", capsize=3, label=f"p={p:g}")
    plt.axhline(0, linewidth=1, linestyle="--")
    plt.xlabel("Resource budget B (arb. units ∝ n·L)")
    plt.ylabel("Control gain ΔI / I [%]")
    plt.title("Monte Carlo gain intervals under finite-shot resampling")
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    FIGS.mkdir(exist_ok=True, parents=True)
    plt.savefig(FIGS / "fig8_mc_gain.png", dpi=300, bbox_inches="tight")
    plt.close()


def mc_coefficients() -> pd.DataFrame:
    coef = pd.read_csv(SUPP / "table_coefficients.csv")
    mu = coef["Estimate"].to_numpy()
    se = coef["Std. Error"].to_numpy()

    rng = np.random.default_rng(SEED + 1)
    N = 200000
    samples = rng.normal(mu, se, size=(N, len(mu)))

    rows = []
    for j, name in enumerate(coef["Parameter"].to_list()):
        s = samples[:, j]
        rows.append(
            dict(
                parameter=name,
                estimate=float(mu[j]),
                std_error=float(se[j]),
                ci_low=float(np.quantile(s, 0.025)),
                ci_high=float(np.quantile(s, 0.975)),
                pr_positive=float(np.mean(s > 0)),
            )
        )

    # Derived tradeoff gamma = a2/a3 (depth vs noise)
    a2 = samples[:, 2]
    a3 = samples[:, 3]
    gamma = a2 / a3
    rows.append(
        dict(
            parameter="gamma=a2/a3",
            estimate=float(np.mean(gamma)),
            std_error=float(np.std(gamma, ddof=1)),
            ci_low=float(np.quantile(gamma, 0.025)),
            ci_high=float(np.quantile(gamma, 0.975)),
            pr_positive=float(np.mean(gamma > 0)),
        )
    )

    out = pd.DataFrame(rows)
    out.to_csv(SUPP / "mc_coefficients_summary.csv", index=False)
    return out


def main() -> None:
    df = mc_gain_intervals()
    plot_gain(df)
    mc_coefficients()
    print("Done. Regenerated Fig. 8 and Monte Carlo summary CSVs.")


if __name__ == "__main__":
    main()



# ==========================
# Extensions: calibration and budget-model robustness
# ==========================

def eta_from_p(p: float) -> float:
    """Depolarizing shrinkage factor for Pauli expectations: eta = 1 - 4p/3."""
    return max(1e-9, 1.0 - 4.0*p/3.0)

def mi_scale_with_eta(mi0: float, eta_hat: float, eta0: float, a3: float) -> float:
    """Local sensitivity model: MI(p) ≈ MI(p0) * (eta/eta0)^{a3}, clipped to [0,1]."""
    if eta0 <= 0:
        return 0.0
    mi = mi0 * (eta_hat/eta0)**a3
    return float(np.clip(mi, 0.0, 0.999999))

def calibration_uncertainty_gain(n_trials: int = 30000, delta_rel: float = 0.10, seed: int = 1234):
    """
    Propagate calibration uncertainty p -> p_hat into the reported baseline/control budget table.
    This is NOT shot noise; it is parameter uncertainty around the nominal p used for each curve.
    """
    rng = np.random.default_rng(seed)

    # Load the budget table that contains baseline and controlled choices + MI values
    df = pd.read_csv("supplement/table_budget_full.csv")

    # Scaling exponent on log(eta) from Table II
    coeff = pd.read_csv("supplement/table_coefficients.csv").set_index("Parameter")["Estimate"].to_dict()
    a3 = float(coeff["log(eta = 1 - 4p/3)"])

    rows = []
    for _, r in df.iterrows():
        p0 = float(r["p"])
        B = float(r["B"])

        # Nominal designs and MI from the budget table (Table III)
        n_b, L_b, mi0_b = parse_baseline(r["Baseline (n,L), I"])
        n_c, L_c, u_c, mi0_c = parse_control(r["Controlled (n,L,u), I"])
        u_b = 0.0

        # Sample calibration-perturbed p_hat (uniform relative band)
        lo = max(0.0, p0*(1.0 - delta_rel))
        hi = p0*(1.0 + delta_rel)
        p_hat = rng.uniform(lo, hi, size=n_trials)

        # Baseline: eta(p_hat)
        eta0_b = eta_from_p(p0*np.exp(-u_b))
        eta_hat_b = 1.0 - 4.0*(p_hat*np.exp(-u_b))/3.0
        eta_hat_b = np.clip(eta_hat_b, 1e-9, 1.0)

        # Control: eta(p_hat e^{-u})
        eta0_c = eta_from_p(p0*np.exp(-u_c))
        eta_hat_c = 1.0 - 4.0*(p_hat*np.exp(-u_c))/3.0
        eta_hat_c = np.clip(eta_hat_c, 1e-9, 1.0)

        mi_b = np.clip(mi0_b * (eta_hat_b/eta0_b)**a3, 0.0, 0.999999)
        mi_c = np.clip(mi0_c * (eta_hat_c/eta0_c)**a3, 0.0, 0.999999)

        gain = (mi_c - mi_b) / np.maximum(mi_b, 1e-9) * 100.0

        rows.append({
            "p": p0, "B": B,
            "gain_median_pct": float(np.median(gain)),
            "gain_p05_pct": float(np.quantile(gain, 0.05)),
            "gain_p95_pct": float(np.quantile(gain, 0.95)),
            "prob_gain_positive": float(np.mean(gain > 0.0)),
        })

    out = pd.DataFrame(rows).sort_values(["p", "B"]).reset_index(drop=True)
    out.to_csv("supplement/p_calibration_gain_summary.csv", index=False)

    # Plot: median gain with 5--95% error bars, grouped by p (similar to Fig. 8 style)
    plt.figure(figsize=(7.0, 4.2))
    for p0 in sorted(out["p"].unique()):
        sub = out[out["p"] == p0].sort_values("B")
        B = sub["B"].values
        med = sub["gain_median_pct"].values
        yerr = np.vstack([med - sub["gain_p05_pct"].values, sub["gain_p95_pct"].values - med])
        plt.errorbar(B, med, yerr=yerr, marker="o", capsize=3, linewidth=1.2, label=f"p={p0:g}")
    plt.xlabel("Resource budget B (arb. units, proportional to nL)")
    plt.ylabel("Control gain in MI [%]")
    plt.title(r"Calibration-uncertainty sensitivity ($p\pm10\%$): control gain")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2, frameon=True)
    plt.tight_layout()
    plt.savefig("figs/fig9_p_calibration_gain.png", dpi=300)
    plt.close()

    return out

def budget_proxy_optimization(seed: int = 123):
    """
    Compare two budget definitions:
      (i) gate-count proxy:  C_gate = n L (1 + c_u u)
      (ii) wall-clock proxy: C_time = n (L + L0) (1 + c_u u)
    using the same scaling-law MI model as in the manuscript.
    """
    coeff = pd.read_csv("supplement/table_coefficients.csv").set_index("Parameter")["Estimate"].to_dict()
    a0 = float(coeff["Intercept"])
    a1 = float(coeff["log(n_qubits)"])
    a2 = float(coeff["log(depth)"])
    a3 = float(coeff["log(eta = 1 - 4p/3)"])

    # Candidate grid for optimization
    n_list = [2, 4, 6, 8]
    L_list = list(range(1, 13))  # allow integer depths 1..12
    u_list = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]

    c_u = 0.6
    L0 = 4  # fixed per-shot overhead in "layer-equivalents" (init/measure/reset proxy)

    def I_pred(n, L, p, u):
        eta = 1.0 - 4.0*(p*np.exp(-u))/3.0
        eta = np.clip(eta, 1e-9, 1.0)
        logI = a0 + a1*np.log(n) + a2*np.log(L) + a3*np.log(eta)
        I = float(np.exp(logI))
        return float(np.clip(I, 0.0, 0.999999))

    def cost_gate(n, L, u):
        return n*L*(1.0 + c_u*u)

    def cost_time(n, L, u):
        return n*(L + L0)*(1.0 + c_u*u)

    p_values = [0.001, 0.005, 0.01, 0.02]
    B_values = [20, 40, 80, 160]

    rows = []
    for p0 in p_values:
        for B in B_values:
            for model_name, cost_fn in [("gate", cost_gate), ("wallclock", cost_time)]:
                # baseline: u=0 only
                best_b = {"I": -1.0}
                for n in n_list:
                    for L in L_list:
                        u = 0.0
                        if cost_fn(n, L, u) <= B + 1e-9:
                            I = I_pred(n, L, p0, u)
                            if I > best_b["I"]:
                                best_b = {"I": I, "n": n, "L": L, "u": u, "cost": cost_fn(n, L, u)}
                # control: allow u>=0
                best_c = {"I": -1.0}
                for n in n_list:
                    for L in L_list:
                        for u in u_list:
                            if cost_fn(n, L, u) <= B + 1e-9:
                                I = I_pred(n, L, p0, u)
                                if I > best_c["I"]:
                                    best_c = {"I": I, "n": n, "L": L, "u": u, "cost": cost_fn(n, L, u)}

                rows.append({
                    "p": p0, "B": B, "budget_model": model_name,
                    "n_baseline": best_b.get("n", np.nan),
                    "L_baseline": best_b.get("L", np.nan),
                    "MI_baseline_pred": best_b.get("I", np.nan),
                    "cost_baseline": best_b.get("cost", np.nan),
                    "n_control": best_c.get("n", np.nan),
                    "L_control": best_c.get("L", np.nan),
                    "u_control": best_c.get("u", np.nan),
                    "MI_control_pred": best_c.get("I", np.nan),
                    "cost_control": best_c.get("cost", np.nan),
                    "gain_pct": (best_c.get("I", np.nan) - best_b.get("I", np.nan)) / max(best_b.get("I", 1e-9), 1e-9) * 100.0
                })

    out = pd.DataFrame(rows).sort_values(["p", "B", "budget_model"]).reset_index(drop=True)
    out.to_csv("supplement/budget_proxy_optimization.csv", index=False)

    # Representative figure: p=0.02 (high noise), compare optimal MI vs B for both budget proxies
    sub = out[(out["p"] == 0.02)].copy()
    sub_gate = sub[sub["budget_model"] == "gate"].sort_values("B")
    sub_time = sub[sub["budget_model"] == "wallclock"].sort_values("B")

    plt.figure(figsize=(7.0, 4.2))
    B = sub_gate["B"].values
    plt.plot(B, sub_gate["MI_baseline_pred"].values, marker="o", linestyle="--", label="Baseline (gate budget)")
    plt.plot(B, sub_gate["MI_control_pred"].values, marker="s", linestyle="-", label="Control (gate budget)")
    plt.plot(B, sub_time["MI_baseline_pred"].values, marker="o", linestyle=":", label="Baseline (wall-clock budget)")
    plt.plot(B, sub_time["MI_control_pred"].values, marker="s", linestyle="-.", label="Control (wall-clock budget)")
    plt.xlabel("Budget B (arb. units)")
    plt.ylabel("Predicted MI [bits]")
    plt.title("Budget-model sensitivity at high noise (p = 0.02)")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8, ncol=2, frameon=True)
    plt.tight_layout()
    plt.savefig("figs/fig10_budget_proxy.png", dpi=300)
    plt.close()

    # Small table for the paper: p=0.02 only
    table_p002 = out[out["p"] == 0.02].copy()
    table_p002.to_csv("supplement/table_budget_proxy_p002.csv", index=False)

    return out

if __name__ == "__main__":
    # Existing outputs (Fig. 8 and Table IV) are generated above.
    # Add robustness extensions used in v3.8+.
    print("\n[Robustness extension] Calibration-uncertainty sensitivity...")
    calibration_uncertainty_gain(n_trials=30000, delta_rel=0.10, seed=1234)

    print("[Robustness extension] Budget-proxy sensitivity...")
    budget_proxy_optimization(seed=123)