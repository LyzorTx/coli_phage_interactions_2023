#!/usr/bin/env python3
"""Case-by-case comparison of two phage-host prediction files.

Reusable diagnostic protocol distilled from SX11/SX12/SX13 analyses. Distinguishes
real-signal / sub-threshold / outlier-cherry-pick reads of small aggregate metric deltas.

Usage:
    python compare_predictions.py <baseline.csv> <candidate.csv> [options]

See SKILL.md for full argument reference.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.metrics import ndcg_score


@dataclass
class Config:
    baseline_path: Path
    candidate_path: Path
    entity_col: str = "bacteria"
    pair_entity_col: str | None = "phage"
    label_col: str = "mlc_score"
    pair_col: str = "pair_id"
    pred_col: str = "predicted_probability"
    focus: list[str] = None
    narrow_rate: float = 0.15
    broad_rate: float = 0.30
    permutations: int = 200
    win_threshold: float = 0.01
    out_path: Path | None = None
    baseline_name: str = "baseline"
    candidate_name: str = "candidate"
    peer_rate_tolerance: float = 0.03


def load_and_merge(cfg: Config) -> pd.DataFrame:
    """Merge two prediction files on pair_col. Each side keeps its predicted_probability."""
    b = pd.read_csv(cfg.baseline_path)
    c = pd.read_csv(cfg.candidate_path)

    for col in (cfg.pair_col, cfg.entity_col, cfg.pred_col):
        if col not in b.columns:
            raise ValueError(f"{cfg.baseline_path}: missing column '{col}'")
        if col not in c.columns:
            raise ValueError(f"{cfg.candidate_path}: missing column '{col}'")

    # Baseline keeps label/entity (plus optional pair_entity_col for named-case tables);
    # candidate only needs pair_id + pred.
    keep_cols = [cfg.pair_col, cfg.entity_col, cfg.label_col, cfg.pred_col]
    if cfg.pair_entity_col and cfg.pair_entity_col in b.columns and cfg.pair_entity_col not in keep_cols:
        keep_cols.append(cfg.pair_entity_col)
    b = b[keep_cols].rename(columns={cfg.pred_col: f"p_{cfg.baseline_name}"})
    c = c[[cfg.pair_col, cfg.pred_col]].rename(columns={cfg.pred_col: f"p_{cfg.candidate_name}"})
    merged = b.merge(c, on=cfg.pair_col, how="inner")

    dropped_b = len(b) - len(merged)
    dropped_c = len(c) - len(merged)
    if dropped_b or dropped_c:
        print(
            f"[warn] Pair-id mismatch: {dropped_b} baseline-only, {dropped_c} candidate-only rows dropped; "
            f"analyzing {len(merged)} shared pairs",
            file=sys.stderr,
        )

    merged[cfg.label_col] = merged[cfg.label_col].fillna(0)
    merged["rel"] = merged[cfg.label_col].astype(int)
    return merged


def per_entity_ndcg(merged: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Per-bacterium nDCG under each prediction file, with lysis rate and positive count."""
    rows = []
    for ent, grp in merged.groupby(cfg.entity_col):
        if grp["rel"].max() == 0 or len(grp) < 2:
            continue
        relevance = grp["rel"].values.reshape(1, -1)
        n_pos = int((grp["rel"] > 0).sum())
        lysis_rate = n_pos / len(grp)
        ndcg_b = ndcg_score(relevance, grp[f"p_{cfg.baseline_name}"].values.reshape(1, -1))
        ndcg_c = ndcg_score(relevance, grp[f"p_{cfg.candidate_name}"].values.reshape(1, -1))
        rows.append(
            {
                cfg.entity_col: ent,
                "n_pos": n_pos,
                "lysis_rate": lysis_rate,
                f"ndcg_{cfg.baseline_name}": ndcg_b,
                f"ndcg_{cfg.candidate_name}": ndcg_c,
                "ndcg_delta": ndcg_c - ndcg_b,
            }
        )
    return pd.DataFrame(rows)


def format_aggregate(df: pd.DataFrame, cfg: Config) -> str:
    out = []
    out.append("## Aggregate per-bacterium nDCG\n")
    b_mean = df[f"ndcg_{cfg.baseline_name}"].mean()
    c_mean = df[f"ndcg_{cfg.candidate_name}"].mean()
    delta_mean = df["ndcg_delta"].mean()
    delta_median = df["ndcg_delta"].median()
    out.append(
        f"- Bacteria analyzed (≥1 positive, ≥2 pairs): **{len(df)}**"
    )
    out.append(f"- Mean nDCG {cfg.baseline_name}: **{b_mean:.4f}**  |  {cfg.candidate_name}: **{c_mean:.4f}**")
    out.append(f"- Mean delta: **{delta_mean:+.4f}** ({delta_mean * 100:+.2f} pp)")
    out.append(f"- Median delta: {delta_median:+.4f} ({delta_median * 100:+.2f} pp)")
    wins = (df["ndcg_delta"] > cfg.win_threshold).sum()
    losses = (df["ndcg_delta"] < -cfg.win_threshold).sum()
    ties = len(df) - wins - losses
    out.append(f"- Wins (Δ ≥ +{cfg.win_threshold:.2f}): **{wins}** | Losses: **{losses}** | Ties: **{ties}**")
    return "\n".join(out) + "\n"


def format_deciles(df: pd.DataFrame, cfg: Config) -> str:
    out = ["## Lysis-rate decile stratification\n"]
    df = df.copy()
    df["decile"] = pd.qcut(df["lysis_rate"], 10, labels=False, duplicates="drop")
    out.append("| Decile | Lysis range | n | Mean Δ nDCG | Median Δ nDCG |")
    out.append("|--------|-------------|---|-------------|----------------|")
    for d in sorted(df["decile"].dropna().unique()):
        sub = df[df["decile"] == d]
        lo, hi = sub["lysis_rate"].min(), sub["lysis_rate"].max()
        mean_d = sub["ndcg_delta"].mean()
        med_d = sub["ndcg_delta"].median()
        out.append(f"| {int(d)} | {lo:.3f}–{hi:.3f} | {len(sub)} | {mean_d:+.4f} | {med_d:+.4f} |")
    return "\n".join(out) + "\n"


def format_narrow_sign_test(df: pd.DataFrame, cfg: Config) -> str:
    narrow = df[df["lysis_rate"] < cfg.narrow_rate]
    broad = df[df["lysis_rate"] >= cfg.broad_rate]
    out = ["## Sign tests\n"]

    def one(label: str, sub: pd.DataFrame) -> str:
        pos = int((sub["ndcg_delta"] > 0).sum())
        n = len(sub)
        if n == 0:
            return f"- **{label}**: 0 bacteria, no test run"
        p = binomtest(pos, n).pvalue
        sig = " 🟢" if p < 0.05 else ""
        return (
            f"- **{label}** (lysis {label.split()[1]}, n={n}): {pos}/{n} positive "
            f"(p = {p:.3f}){sig}  |  mean Δ = {sub['ndcg_delta'].mean():+.4f}"
        )

    out.append(one(f"Narrow <{cfg.narrow_rate:.0%}", narrow))
    out.append(one(f"Broad  ≥{cfg.broad_rate:.0%}", broad))
    return "\n".join(out) + "\n"


def permutation_test(merged: pd.DataFrame, cfg: Config, n_perms: int) -> tuple[float, float, float, float]:
    """Swap each pair's baseline/candidate preds 50/50 and recompute aggregate mean Δ nDCG.

    Returns (actual_delta, null_mean, null_std, frac_as_extreme).
    """
    rng = np.random.default_rng(42)
    actual_df = per_entity_ndcg(merged, cfg)
    actual_delta = float(actual_df["ndcg_delta"].mean())

    base_pred = merged[f"p_{cfg.baseline_name}"].to_numpy(copy=True)
    cand_pred = merged[f"p_{cfg.candidate_name}"].to_numpy(copy=True)
    n = len(merged)

    null_deltas = np.empty(n_perms, dtype=float)
    for i in range(n_perms):
        swap = rng.random(n) > 0.5
        perm = merged.copy()
        new_base = np.where(swap, cand_pred, base_pred)
        new_cand = np.where(swap, base_pred, cand_pred)
        perm[f"p_{cfg.baseline_name}"] = new_base
        perm[f"p_{cfg.candidate_name}"] = new_cand
        perm_df = per_entity_ndcg(perm, cfg)
        null_deltas[i] = perm_df["ndcg_delta"].mean()

    null_mean = float(null_deltas.mean())
    null_std = float(null_deltas.std())
    frac = float((np.abs(null_deltas) >= abs(actual_delta)).mean())
    return actual_delta, null_mean, null_std, frac


def format_permutation(result: tuple[float, float, float, float], n_perms: int) -> str:
    actual, m, s, frac = result
    out = ["## Permutation test (random pair-level swap)\n"]
    out.append(f"- Actual aggregate mean Δ nDCG: **{actual:+.4f}**")
    out.append(f"- Null distribution (n={n_perms}): mean {m:+.4f}, std {s:.4f}")
    out.append(f"- Fraction of perms as extreme as observed: **{frac * 100:.1f}%**")
    if frac >= 0.1:
        out.append("- 🔴 Signal **indistinguishable from noise** (≥10% of random swaps are as extreme)")
    elif frac >= 0.05:
        out.append("- 🟡 Marginal (5–10% of swaps as extreme) — directional but not significant")
    else:
        out.append("- 🟢 Signal **distinguishable from noise** (<5% of swaps as extreme)")
    return "\n".join(out) + "\n"


def format_focus_case(merged: pd.DataFrame, per_ent: pd.DataFrame, cfg: Config, bact: str) -> str:
    sub = per_ent[per_ent[cfg.entity_col] == bact]
    if sub.empty:
        return f"\n### {bact}: not in per-entity table (no positives or <2 pairs)\n"
    row = sub.iloc[0]
    out = [f"\n### Spotlight: {bact}\n"]
    out.append(f"- n_pos = {int(row['n_pos'])}, lysis_rate = {row['lysis_rate']:.3f}")
    out.append(
        f"- nDCG {cfg.baseline_name}: {row[f'ndcg_{cfg.baseline_name}']:.4f}  →  "
        f"{cfg.candidate_name}: {row[f'ndcg_{cfg.candidate_name}']:.4f}  "
        f"(Δ {row['ndcg_delta']:+.4f})"
    )

    # Peer comparison: bacteria within ±peer_rate_tolerance of this bacterium's lysis rate
    rate = row["lysis_rate"]
    peers = per_ent[
        (per_ent["lysis_rate"] >= rate - cfg.peer_rate_tolerance)
        & (per_ent["lysis_rate"] <= rate + cfg.peer_rate_tolerance)
    ]
    if len(peers) > 1:
        peer_mean = peers["ndcg_delta"].mean()
        peer_median = peers["ndcg_delta"].median()
        percentile = float((peers["ndcg_delta"] < row["ndcg_delta"]).mean() * 100)
        out.append(
            f"- Peer comparison (±{cfg.peer_rate_tolerance:.2f} lysis rate, n={len(peers)}): "
            f"peer mean Δ = {peer_mean:+.4f}, median = {peer_median:+.4f}, this bacterium at **{percentile:.0f}th percentile**"
        )
        if percentile > 85 or percentile < 15:
            out.append(
                "- ⚠️  **Outlier relative to peers** — do not generalize this bacterium's Δ to narrow-host class"
            )
        else:
            out.append("- Typical for its peer group")

    # Per-pair rank changes for positives
    pairs = merged[merged[cfg.entity_col] == bact].copy()
    pos = pairs[pairs["rel"] > 0].sort_values(f"p_{cfg.baseline_name}", ascending=False)
    if len(pos) > 0:
        pair_col_name = cfg.pair_entity_col if cfg.pair_entity_col and cfg.pair_entity_col in pos.columns else None
        header = "Phage" if pair_col_name == "phage" else (pair_col_name or "Pair entity")
        out.append(f"\n| {header} | MLC | Baseline rank | Candidate rank | Δ rank |")
        out.append("|-------|-----|---------------|----------------|--------|")
        pred_b = pairs[f"p_{cfg.baseline_name}"].to_numpy()
        pred_c = pairs[f"p_{cfg.candidate_name}"].to_numpy()
        for _, r in pos.iterrows():
            rb = int((pred_b >= r[f"p_{cfg.baseline_name}"]).sum())
            rc = int((pred_c >= r[f"p_{cfg.candidate_name}"]).sum())
            name = r[pair_col_name] if pair_col_name else "?"
            out.append(f"| {name} | {int(r['rel'])} | {rb} | {rc} | {rc - rb:+d} |")
    return "\n".join(out) + "\n"


def format_top3(merged: pd.DataFrame, cfg: Config) -> str:
    out = ["## Top-3 hit rate (any positive in top-3 predictions)\n"]
    hits_b, hits_c, n = 0, 0, 0
    for ent, grp in merged.groupby(cfg.entity_col):
        if grp["rel"].max() == 0:
            continue
        n += 1
        top3_b = grp.nlargest(3, f"p_{cfg.baseline_name}")
        top3_c = grp.nlargest(3, f"p_{cfg.candidate_name}")
        if (top3_b["rel"] > 0).any():
            hits_b += 1
        if (top3_c["rel"] > 0).any():
            hits_c += 1
    out.append(f"- {cfg.baseline_name}: **{hits_b}/{n} = {100 * hits_b / n:.1f}%**")
    out.append(f"- {cfg.candidate_name}: **{hits_c}/{n} = {100 * hits_c / n:.1f}%**")
    out.append(f"- Net: **{hits_c - hits_b:+d} bacteria** (positive = rescued, negative = lost)")
    return "\n".join(out) + "\n"


def build_decision_headline(
    per_ent: pd.DataFrame,
    perm_result: tuple[float, float, float, float],
    narrow_df: pd.DataFrame,
    cfg: Config,
    focus_bacts: Iterable[str],
) -> str:
    actual, _, _, frac = perm_result
    delta_pp = actual * 100
    out = ["# Headline\n"]
    out.append(f"- Aggregate mean Δ nDCG: **{delta_pp:+.2f} pp**")
    if frac < 0.05:
        sig = "distinguishable from noise"
    elif frac < 0.10:
        sig = "marginal (5–10% permutation p)"
    else:
        sig = "indistinguishable from noise"
    out.append(f"- Permutation significance: **{sig}** (p ≈ {frac:.2f})")
    if len(narrow_df) > 0:
        pos = int((narrow_df["ndcg_delta"] > 0).sum())
        n = len(narrow_df)
        p = binomtest(pos, n).pvalue
        out.append(f"- Narrow hosts (lysis <{cfg.narrow_rate:.0%}): {pos}/{n} positive, sign test p = {p:.3f}")
    out.append("\n# Decision\n")
    if frac < 0.05 and abs(delta_pp) >= 2:
        out.append("- **Real signal.** Clears conventional +2 pp bar AND survives permutation test.")
    elif frac < 0.05:
        out.append(
            "- **Sub-threshold but real signal.** Survives permutation test but below +2 pp. "
            "Note in notebook; do not adopt as production change."
        )
    elif frac < 0.10:
        out.append(
            "- **Marginal.** Directionally coherent but statistically weak. "
            "Record as observation; consider as hypothesis for future work."
        )
    else:
        out.append(
            "- **Null.** Aggregate delta indistinguishable from random prediction-swap noise. "
            "Record as validated null."
        )
    for b in focus_bacts or []:
        row = per_ent[per_ent[cfg.entity_col] == b]
        if row.empty:
            continue
        delta = float(row["ndcg_delta"].iloc[0])
        if abs(delta) > 0.02:
            out.append(f"- Spotlight **{b}** Δ {delta:+.4f}: check peer comparison before citing as rescue / regression.")
    return "\n".join(out) + "\n"


def parse_args(argv: list[str] | None = None) -> Config:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("baseline", type=Path, help="Baseline predictions CSV")
    p.add_argument("candidate", type=Path, help="Candidate predictions CSV")
    p.add_argument("--entity-col", default="bacteria")
    p.add_argument(
        "--pair-entity-col",
        default="phage",
        help="Name of the 'other side' of each interaction (used in spotlight per-pair table). Default: phage. "
        "Set to empty string to disable.",
    )
    p.add_argument("--label-col", default="mlc_score", help="Graded relevance column (MLC 0-3 by default)")
    p.add_argument("--pair-col", default="pair_id")
    p.add_argument("--pred-col", default="predicted_probability")
    p.add_argument("--focus", default="", help="Comma-separated list of entity names to spotlight")
    p.add_argument("--narrow-rate", type=float, default=0.15)
    p.add_argument("--broad-rate", type=float, default=0.30)
    p.add_argument("--permutations", type=int, default=200)
    p.add_argument("--win-threshold", type=float, default=0.01)
    p.add_argument("--baseline-name", default="baseline")
    p.add_argument("--candidate-name", default="candidate")
    p.add_argument("--peer-rate-tolerance", type=float, default=0.03)
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args(argv)
    return Config(
        baseline_path=a.baseline,
        candidate_path=a.candidate,
        entity_col=a.entity_col,
        pair_entity_col=a.pair_entity_col or None,
        label_col=a.label_col,
        pair_col=a.pair_col,
        pred_col=a.pred_col,
        focus=[x.strip() for x in a.focus.split(",") if x.strip()],
        narrow_rate=a.narrow_rate,
        broad_rate=a.broad_rate,
        permutations=a.permutations,
        win_threshold=a.win_threshold,
        baseline_name=a.baseline_name,
        candidate_name=a.candidate_name,
        peer_rate_tolerance=a.peer_rate_tolerance,
        out_path=a.out,
    )


def main(argv: list[str] | None = None) -> None:
    cfg = parse_args(argv)
    merged = load_and_merge(cfg)
    per_ent = per_entity_ndcg(merged, cfg)

    narrow_df = per_ent[per_ent["lysis_rate"] < cfg.narrow_rate]

    sections: list[str] = []
    sections.append(f"# Case-by-case comparison: {cfg.baseline_name} vs {cfg.candidate_name}\n")
    sections.append(f"Baseline file: `{cfg.baseline_path}`")
    sections.append(f"Candidate file: `{cfg.candidate_path}`\n")

    sections.append(format_aggregate(per_ent, cfg))
    sections.append(format_deciles(per_ent, cfg))
    sections.append(format_narrow_sign_test(per_ent, cfg))
    sections.append(format_top3(merged, cfg))

    perm_result = permutation_test(merged, cfg, cfg.permutations)
    sections.append(format_permutation(perm_result, cfg.permutations))

    if cfg.focus:
        sections.append("## Named-case spotlight\n")
        for b in cfg.focus:
            sections.append(format_focus_case(merged, per_ent, cfg, b))

    sections.append(build_decision_headline(per_ent, perm_result, narrow_df, cfg, cfg.focus))

    report = "\n".join(sections)
    print(report)
    if cfg.out_path:
        cfg.out_path.parent.mkdir(parents=True, exist_ok=True)
        cfg.out_path.write_text(report, encoding="utf-8")
        print(f"\n[wrote] {cfg.out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
