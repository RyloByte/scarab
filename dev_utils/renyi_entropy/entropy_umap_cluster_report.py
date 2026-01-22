#!/usr/bin/env python3
"""
Path: scripts/entropy_umap_cluster_report.py

What this does:
- Reads your table (TSV/CSV) containing sample_id, sample_type, cluster, Renyi_Entropy, u0, u1 (+ optionally alpha)
- Produces:
  - Agreement stats between sample_type and cluster (ARI/NMI/chi-square + contingency tables)
  - UMAP-space separation stats (silhouette + centroid permutation test)
  - Group differences in Renyi_Entropy (Kruskal + pairwise Mann–Whitney + FDR)
  - **NEW**: PERMANOVA on the *entropy feature matrix* (wide: sample_id x alpha if alpha exists; otherwise scalar Renyi_Entropy)

Run:
  python scripts/entropy_umap_cluster_report.py \
    --input data/table.tsv \
    --outdir results/entropy_umap_report \
    --sep "\t" \
    --n-perm 5000 \
    --seed 0
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from statsmodels.stats.multitest import multipletests


# -----------------------------
# IO + cleaning
# -----------------------------

def _read_table(path: str, sep: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input not found: {p}")
    if sep.lower() == "auto":
        sep_use = "\t" if p.suffix.lower() in [".tsv", ".tab"] else ","
    else:
        sep_use = sep
    return pd.read_csv(p, sep=sep_use)


def _ensure_cols(df: pd.DataFrame, required: list[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _clean_base(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # Required columns
    df["sample_id"] = df["sample_id"].astype(str)
    df["sample_type"] = df["sample_type"].astype(str)
    df["cluster"] = df["cluster"].astype(str)

    for c in ["Renyi_Entropy", "u0", "u1"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # alpha is optional (if present, used for PERMANOVA feature matrix)
    if "alpha" in df.columns:
        df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")

    # Drop missing essentials for the core report (UMAP + Renyi tests)
    df = df.dropna(subset=["sample_id", "sample_type", "cluster", "Renyi_Entropy", "u0", "u1"])
    return df


def _reduce_to_one_row_per_sample(df: pd.DataFrame) -> pd.DataFrame:
    """
    If the input is long-form (multiple rows per sample_id, e.g. multiple alpha),
    this collapses to one row per sample_id for the UMAP-based plots/stats.

    Strategy:
    - take first u0/u1/sample_type/cluster seen for each sample_id
    - take median Renyi_Entropy across rows for that sample_id (robust)
    """
    cols_take_first = ["sample_type", "cluster", "u0", "u1"]
    agg = {c: "first" for c in cols_take_first}
    agg["Renyi_Entropy"] = "median"

    out = (
        df.groupby("sample_id", as_index=False)
          .agg(agg)
    )
    return out


# -----------------------------
# Pairwise nonparametrics
# -----------------------------

def _pairwise_mannwhitney(values: pd.Series, groups: pd.Series) -> pd.DataFrame:
    g = groups.astype(str)
    cats = sorted(g.unique())
    rows = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            a = values[g == cats[i]].values
            b = values[g == cats[j]].values
            if len(a) < 2 or len(b) < 2:
                u = np.nan
                p = np.nan
            else:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rows.append(
                {
                    "group_a": cats[i],
                    "group_b": cats[j],
                    "U": u,
                    "p_raw": p,
                    "n_a": len(a),
                    "n_b": len(b),
                    "median_a": float(np.nanmedian(a)) if len(a) else np.nan,
                    "median_b": float(np.nanmedian(b)) if len(b) else np.nan,
                }
            )
    out = pd.DataFrame(rows)
    if out["p_raw"].notna().any():
        mask = out["p_raw"].notna()
        rej, p_adj, _, _ = multipletests(out.loc[mask, "p_raw"].values, method="fdr_bh")
        out.loc[mask, "p_fdr_bh"] = p_adj
        out.loc[mask, "reject_fdr_0.05"] = rej
    else:
        out["p_fdr_bh"] = np.nan
        out["reject_fdr_0.05"] = np.nan
    return out.sort_values(["p_fdr_bh", "p_raw"], na_position="last")


# -----------------------------
# UMAP-space separation tests
# -----------------------------

def _centroid_separation_permtest(X: np.ndarray, labels: np.ndarray, n_perm: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)

    def metric(lab: np.ndarray) -> float:
        labs = pd.Series(lab).astype(str).values
        uniq = np.unique(labs)
        if len(uniq) < 2:
            return np.nan
        cents = []
        for u in uniq:
            pts = X[labs == u]
            if pts.shape[0] == 0:
                continue
            cents.append(pts.mean(axis=0))
        cents = np.vstack(cents)
        d = []
        for i in range(cents.shape[0]):
            for j in range(i + 1, cents.shape[0]):
                d.append(np.linalg.norm(cents[i] - cents[j]))
        return float(np.mean(d)) if d else np.nan

    obs = metric(labels)
    null = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        null[k] = metric(rng.permutation(labels))

    p = (np.sum(null >= obs) + 1.0) / (n_perm + 1.0) if np.isfinite(obs) else np.nan
    return {
        "obs_mean_centroid_dist": obs,
        "p_perm": p,
        "n_perm": n_perm,
        "null_mean": float(np.nanmean(null)),
        "null_sd": float(np.nanstd(null)),
    }


def _savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


# -----------------------------
# PERMANOVA (no extra deps)
# -----------------------------

def _euclidean_distance_matrix(X: np.ndarray) -> np.ndarray:
    """
    Pairwise Euclidean distances for an (n x p) matrix.
    Returns (n x n) distance matrix.
    """
    # (x - y)^2 = x^2 + y^2 - 2xy
    G = X @ X.T
    sq = np.diag(G).reshape(-1, 1)
    D2 = np.maximum(sq + sq.T - 2.0 * G, 0.0)
    return np.sqrt(D2)


def _permanova(D: np.ndarray, groups: np.ndarray, n_perm: int, seed: int) -> dict:
    """
    One-way PERMANOVA on a distance matrix D with grouping labels.

    Implements the common Anderson-style pseudo-F using distance-based sums of squares:
      SSt = sum_{i<j} d_ij^2 / n
      SSw = sum_g sum_{i<j in g} d_ij^2 / n_g
      SSb = SSt - SSw
      F = (SSb/(k-1)) / (SSw/(n-k))

    Permutation test by shuffling group labels.

    Returns: pseudo-F, p-value, dfs, SS terms, n/k
    """
    rng = np.random.default_rng(seed)
    labs = pd.Series(groups).astype(str).values
    n = D.shape[0]
    uniq = np.unique(labs)
    k = len(uniq)

    if k < 2 or n <= k:
        return {
            "n": n, "k": k,
            "pseudo_F": np.nan,
            "p_perm": np.nan,
            "df_between": np.nan,
            "df_within": np.nan,
            "SS_total": np.nan,
            "SS_within": np.nan,
            "SS_between": np.nan,
            "n_perm": n_perm,
        }

    D2 = D ** 2

    def ss_total() -> float:
        return float(np.sum(np.triu(D2, 1)) / n)

    def ss_within(l: np.ndarray) -> float:
        out = 0.0
        for u in np.unique(l):
            idx = np.where(l == u)[0]
            ng = len(idx)
            if ng < 2:
                continue
            out += np.sum(np.triu(D2[np.ix_(idx, idx)], 1)) / ng
        return float(out)

    SSt = ss_total()
    SSw = ss_within(labs)
    SSb = SSt - SSw

    dfb = k - 1
    dfw = n - k
    MSb = SSb / dfb
    MSw = SSw / dfw
    F_obs = MSb / MSw if MSw > 0 else np.nan

    # permutation null
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        perm = rng.permutation(labs)
        SSw_p = ss_within(perm)
        SSb_p = SSt - SSw_p
        MSb_p = SSb_p / dfb
        MSw_p = SSw_p / dfw
        null[i] = MSb_p / MSw_p if MSw_p > 0 else np.nan

    p = (np.sum(null >= F_obs) + 1.0) / (n_perm + 1.0) if np.isfinite(F_obs) else np.nan

    return {
        "n": n, "k": k,
        "pseudo_F": float(F_obs),
        "p_perm": float(p),
        "df_between": dfb,
        "df_within": dfw,
        "SS_total": float(SSt),
        "SS_within": float(SSw),
        "SS_between": float(SSb),
        "n_perm": n_perm,
        "null_mean": float(np.nanmean(null)),
        "null_sd": float(np.nanstd(null)),
    }


def _build_entropy_feature_matrix(df_raw: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Builds an entropy feature matrix for PERMANOVA.

    If 'alpha' exists and there are multiple rows per sample_id:
      X = pivot(sample_id x alpha) with Renyi_Entropy values
    Else:
      X = single-column matrix with Renyi_Entropy median per sample_id

    Returns:
      X_df: index=sample_id, columns=features
      meta: sample_id, sample_type, cluster, u0, u1 (1 row per sample_id)
    """
    # meta: first seen per sample_id (labels/coords)
    meta = (
        df_raw.sort_values(["sample_id"])
              .groupby("sample_id", as_index=False)
              .agg({
                  "sample_type": "first",
                  "cluster": "first",
                  "u0": "first",
                  "u1": "first",
              })
    )

    if "alpha" in df_raw.columns and df_raw["alpha"].notna().any():
        # wide matrix: sample_id x alpha
        X = (
            df_raw.dropna(subset=["alpha", "Renyi_Entropy"])
                 .pivot_table(index="sample_id", columns="alpha", values="Renyi_Entropy", aggfunc="median")
                 .sort_index(axis=1)
        )
        # rename columns to stable strings
        X.columns = [f"alpha_{a:g}" for a in X.columns]
    else:
        # scalar entropy per sample
        X = (
            df_raw.groupby("sample_id")["Renyi_Entropy"]
                  .median()
                  .to_frame(name="Renyi_Entropy")
        )

    # Align meta and X; keep intersection only
    common = sorted(set(meta["sample_id"]).intersection(set(X.index)))
    meta = meta[meta["sample_id"].isin(common)].set_index("sample_id").loc[common].reset_index()
    X = X.loc[common]

    # Impute missing features with column medians (PERMANOVA requires full matrix)
    X_imp = X.copy()
    for c in X_imp.columns:
        med = float(np.nanmedian(X_imp[c].values))
        X_imp[c] = X_imp[c].fillna(med)

    return X_imp, meta


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Input TSV/CSV with required columns.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--sep", default="auto", help=r"Delimiter: '\t', ',', or 'auto'.")
    ap.add_argument("--n-perm", type=int, default=5000, help="Permutations for permutation-based tests.")
    ap.add_argument("--seed", type=int, default=0, help="Random seed.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    _mkdir(outdir)
    _mkdir(outdir / "tables")
    _mkdir(outdir / "plots")

    # Read + clean
    df_raw = _read_table(args.input, args.sep)
    _ensure_cols(df_raw, ["sample_id", "sample_type", "cluster", "Renyi_Entropy", "u0", "u1"])
    df_raw = _clean_base(df_raw)

    # Build per-sample table for the original plots/tests
    df = _reduce_to_one_row_per_sample(df_raw)

    # -----------------------------
    # 1) Agreement between labelings
    # -----------------------------
    ct = pd.crosstab(df["sample_type"], df["cluster"])
    ct.to_csv(outdir / "tables" / "contingency_sample_type_vs_cluster.csv")

    ari = adjusted_rand_score(df["sample_type"].astype(str), df["cluster"].astype(str))
    nmi = normalized_mutual_info_score(df["sample_type"].astype(str), df["cluster"].astype(str))

    chi2, chi2_p, dof, expected = stats.chi2_contingency(ct.values)
    expected_df = pd.DataFrame(expected, index=ct.index, columns=ct.columns)
    expected_df.to_csv(outdir / "tables" / "contingency_expected_counts.csv")

    agreement_summary = pd.DataFrame([{
        "n_samples": df.shape[0],
        "n_sample_types": df["sample_type"].nunique(),
        "n_clusters": df["cluster"].nunique(),
        "ARI(sample_type,cluster)": ari,
        "NMI(sample_type,cluster)": nmi,
        "chi2_p_value": chi2_p,
        "chi2_stat": chi2,
        "chi2_dof": dof,
    }])
    agreement_summary.to_csv(outdir / "tables" / "agreement_summary.csv", index=False)

    plt.figure(figsize=(8, 6))
    plt.imshow(ct.values, aspect="auto")
    plt.xticks(np.arange(ct.shape[1]), ct.columns, rotation=45, ha="right")
    plt.yticks(np.arange(ct.shape[0]), ct.index)
    plt.colorbar(label="count")
    plt.title("sample_type vs cluster (counts)")
    _savefig(outdir / "plots" / "heatmap_sample_type_vs_cluster.png")

    # -----------------------------
    # 2) Separation in UMAP space
    # -----------------------------
    X_umap = df[["u0", "u1"]].to_numpy()

    def safe_silhouette(labels: pd.Series) -> float:
        labs = labels.astype(str).values
        if len(np.unique(labs)) < 2:
            return np.nan
        return float(silhouette_score(X_umap, labs, metric="euclidean"))

    sil_type = safe_silhouette(df["sample_type"])
    sil_clust = safe_silhouette(df["cluster"])

    sep_type = _centroid_separation_permtest(X_umap, df["sample_type"].astype(str).values, args.n_perm, args.seed)
    sep_clust = _centroid_separation_permtest(X_umap, df["cluster"].astype(str).values, args.n_perm, args.seed)

    separation_summary = pd.DataFrame([{
        "silhouette_sample_type": sil_type,
        "silhouette_cluster": sil_clust,
        **{f"sample_type_{k}": v for k, v in sep_type.items()},
        **{f"cluster_{k}": v for k, v in sep_clust.items()},
    }])
    separation_summary.to_csv(outdir / "tables" / "embedding_separation_summary.csv", index=False)

    plt.figure(figsize=(7, 6))
    for lab in sorted(df["sample_type"].unique(), key=str):
        sub = df[df["sample_type"] == lab]
        plt.scatter(sub["u0"], sub["u1"], s=18, alpha=0.7, label=str(lab))
    plt.xlabel("u0")
    plt.ylabel("u1")
    plt.title("UMAP embedding colored by sample_type")
    plt.legend(markerscale=1.2, bbox_to_anchor=(1.02, 1), loc="upper left")
    _savefig(outdir / "plots" / "umap_by_sample_type.png")

    plt.figure(figsize=(7, 6))
    for lab in sorted(df["cluster"].unique(), key=str):
        sub = df[df["cluster"] == lab]
        plt.scatter(sub["u0"], sub["u1"], s=18, alpha=0.7, label=str(lab))
    plt.xlabel("u0")
    plt.ylabel("u1")
    plt.title("UMAP embedding colored by cluster")
    plt.legend(markerscale=1.2, bbox_to_anchor=(1.02, 1), loc="upper left")
    _savefig(outdir / "plots" / "umap_by_cluster.png")

    # -----------------------------
    # 3) Renyi_Entropy differences by group
    # -----------------------------
    def kruskal_by(label_col: str) -> dict:
        groups = []
        for lab in sorted(df[label_col].unique(), key=str):
            vals = df.loc[df[label_col] == lab, "Renyi_Entropy"].values
            if len(vals) > 0:
                groups.append(vals)
        if len(groups) < 2:
            return {"kruskal_H": np.nan, "kruskal_p": np.nan}
        H, p = stats.kruskal(*groups)
        return {"kruskal_H": float(H), "kruskal_p": float(p)}

    kw_type = kruskal_by("sample_type")
    kw_clust = kruskal_by("cluster")

    renyi_summary = pd.DataFrame([{
        "Renyi_vs_sample_type_kruskal_H": kw_type["kruskal_H"],
        "Renyi_vs_sample_type_kruskal_p": kw_type["kruskal_p"],
        "Renyi_vs_cluster_kruskal_H": kw_clust["kruskal_H"],
        "Renyi_vs_cluster_kruskal_p": kw_clust["kruskal_p"],
    }])
    renyi_summary.to_csv(outdir / "tables" / "renyi_group_tests.csv", index=False)

    pw_type = _pairwise_mannwhitney(df["Renyi_Entropy"], df["sample_type"])
    pw_clust = _pairwise_mannwhitney(df["Renyi_Entropy"], df["cluster"])
    pw_type.to_csv(outdir / "tables" / "pairwise_Renyi_by_sample_type_mannwhitney_fdr.csv", index=False)
    pw_clust.to_csv(outdir / "tables" / "pairwise_Renyi_by_cluster_mannwhitney_fdr.csv", index=False)

    plt.figure(figsize=(9, 5))
    order = sorted(df["sample_type"].unique(), key=str)
    data = [df.loc[df["sample_type"] == lab, "Renyi_Entropy"].values for lab in order]
    plt.boxplot(data, labels=order, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Renyi_Entropy")
    plt.title("Renyi_Entropy by sample_type")
    _savefig(outdir / "plots" / "renyi_by_sample_type_boxplot.png")

    plt.figure(figsize=(9, 5))
    order = sorted(df["cluster"].unique(), key=str)
    data = [df.loc[df["cluster"] == lab, "Renyi_Entropy"].values for lab in order]
    plt.boxplot(data, labels=order, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Renyi_Entropy")
    plt.title("Renyi_Entropy by cluster")
    _savefig(outdir / "plots" / "renyi_by_cluster_boxplot.png")

    # -----------------------------
    # 5) PERMANOVA on entropy feature matrix (reviewer-friendly)
    #    Ctrl+F: "5) PERMANOVA"
    # -----------------------------
    X_entropy, meta = _build_entropy_feature_matrix(df_raw)
    D = _euclidean_distance_matrix(X_entropy.to_numpy())

    perm_type = _permanova(D, meta["sample_type"].astype(str).values, args.n_perm, args.seed)
    perm_clust = _permanova(D, meta["cluster"].astype(str).values, args.n_perm, args.seed)

    permanova_summary = pd.DataFrame([{
        "feature_matrix_shape": f"{X_entropy.shape[0]}x{X_entropy.shape[1]}",
        **{f"PERMANOVA_sample_type_{k}": v for k, v in perm_type.items()},
        **{f"PERMANOVA_cluster_{k}": v for k, v in perm_clust.items()},
    }])
    permanova_summary.to_csv(outdir / "tables" / "permanova_entropy_features.csv", index=False)

    # -----------------------------
    # SUMMARY.txt
    # -----------------------------
    summary_txt = []
    summary_txt.append(f"n_samples = {df.shape[0]}")
    summary_txt.append(f"ARI(sample_type, cluster) = {ari:.4f}")
    summary_txt.append(f"NMI(sample_type, cluster) = {nmi:.4f}")
    summary_txt.append(f"Chi-square p(sample_type vs cluster) = {chi2_p:.3e}")
    summary_txt.append(f"Silhouette(sample_type) = {sil_type:.4f}" if np.isfinite(sil_type) else "Silhouette(sample_type) = NA")
    summary_txt.append(f"Silhouette(cluster) = {sil_clust:.4f}" if np.isfinite(sil_clust) else "Silhouette(cluster) = NA")
    summary_txt.append(f"Centroid sep perm p(sample_type) = {sep_type['p_perm']:.3e}" if np.isfinite(sep_type["p_perm"]) else "Centroid sep perm p(sample_type) = NA")
    summary_txt.append(f"Centroid sep perm p(cluster) = {sep_clust['p_perm']:.3e}" if np.isfinite(sep_clust["p_perm"]) else "Centroid sep perm p(cluster) = NA")
    summary_txt.append(f"Kruskal p(Renyi ~ sample_type) = {kw_type['kruskal_p']:.3e}" if np.isfinite(kw_type["kruskal_p"]) else "Kruskal p(Renyi ~ sample_type) = NA")
    summary_txt.append(f"Kruskal p(Renyi ~ cluster) = {kw_clust['kruskal_p']:.3e}" if np.isfinite(kw_clust["kruskal_p"]) else "Kruskal p(Renyi ~ cluster) = NA")

    summary_txt.append(f"PERMANOVA pseudo-F (sample_type) = {perm_type['pseudo_F']:.4f}" if np.isfinite(perm_type["pseudo_F"]) else "PERMANOVA pseudo-F (sample_type) = NA")
    summary_txt.append(f"PERMANOVA p (sample_type) = {perm_type['p_perm']:.3e}" if np.isfinite(perm_type["p_perm"]) else "PERMANOVA p (sample_type) = NA")
    summary_txt.append(f"PERMANOVA pseudo-F (cluster) = {perm_clust['pseudo_F']:.4f}" if np.isfinite(perm_clust["pseudo_F"]) else "PERMANOVA pseudo-F (cluster) = NA")
    summary_txt.append(f"PERMANOVA p (cluster) = {perm_clust['p_perm']:.3e}" if np.isfinite(perm_clust["p_perm"]) else "PERMANOVA p (cluster) = NA")

    (outdir / "SUMMARY.txt").write_text("\n".join(summary_txt) + "\n")

    print(f"[OK] Wrote report to: {outdir}")
    print(f"[OK] Key summary: {outdir / 'SUMMARY.txt'}")
    print(f"[OK] PERMANOVA table: {outdir / 'tables' / 'permanova_entropy_features.csv'}")


if __name__ == "__main__":
    main()
