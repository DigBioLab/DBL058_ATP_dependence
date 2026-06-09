#!/usr/bin/env python3
"""
Plot ATP conformer cluster embedding and bar chart from saved analysis results.
Requires: numpy, scipy, matplotlib, rdkit

Usage:
    python scripts/plot_clusters.py [--output_dir output] [--plot_dir pictures] [--cutoff 1.5]
"""

import os, json, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from rdkit.ML.Cluster import Butina
from scipy.sparse.linalg import LinearOperator, eigsh

SMALL_THRESHOLD = 10


def mds_embedding(rmsd_mat, cache_path):
    if Path(cache_path).exists():
        print(f"  Loading cached MDS embedding from {cache_path}")
        return np.load(cache_path)

    print("  Computing classical MDS (PCoA) via Lanczos eigensolver …")
    print("  (this takes a few minutes for ~8k conformers)")
    D2 = rmsd_mat.astype(np.float32) ** 2
    n = D2.shape[0]

    # B = -½ H D² H   via matrix-vector products (avoids forming B explicitly)
    def matvec(v):
        w = v - v.mean()
        Bv = D2 @ w.astype(np.float32)
        Bv -= Bv.mean()
        return -0.5 * Bv.astype(np.float64)

    B_op = LinearOperator((n, n), matvec=matvec, dtype=np.float64)
    vals, vecs = eigsh(B_op, k=2, which="LM", tol=1e-4, maxiter=500)

    # Sort descending
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    vals = np.maximum(vals, 0.0)
    embedding = vecs * np.sqrt(vals)
    np.save(cache_path, embedding)
    print(f"  Saved to {cache_path}")
    return embedding


def plot_embedding(embedding, labels, clusters, big_indices, output_path):
    cmap = plt.get_cmap("tab10" if len(big_indices) <= 10 else "tab20")
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor("#f9f9f9")
    ax.set_facecolor("#f9f9f9")

    # Small / unassigned: grey background
    small_mask = ~np.isin(labels, big_indices)
    ax.scatter(
        embedding[small_mask, 0], embedding[small_mask, 1],
        c="#b0b0b0", s=4, alpha=0.25, linewidths=0,
        label="Small clusters", rasterized=True, zorder=1,
    )

    # Big clusters: coloured, front
    for rank_i, ci in enumerate(big_indices):
        mask = labels == ci
        color = cmap(rank_i / max(len(big_indices) - 1, 1))
        ax.scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=[color], s=10, alpha=0.75, linewidths=0,
            label=f"Cluster {rank_i + 1}  (n={len(clusters[ci]):,})",
            rasterized=True, zorder=2 + rank_i,
        )

    n_total = len(labels)
    n_clusters = len(clusters)
    cutoff = 1.5  # stored in clusters.json but passed via caller
    ax.set_title(
        f"ATP Conformer Cluster Embedding — Classical MDS on RMSD\n"
        f"{n_total:,} conformers  ·  {n_clusters} clusters",
        fontsize=12,
    )
    ax.set_xlabel("MDS 1", fontsize=10)
    ax.set_ylabel("MDS 2", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, markerscale=2.5,
              framealpha=0.9, edgecolor="#cccccc")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False,
                   labelleft=False, labelbottom=False)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {output_path}")


def plot_bar(clusters, output_path):
    sizes = [len(cl) for cl in clusters]
    big_bars  = [(i + 1, s) for i, s in enumerate(sizes) if s >= SMALL_THRESHOLD]
    small_total = sum(s for s in sizes if s < SMALL_THRESHOLD)
    n_small_clusters = sum(1 for s in sizes if s < SMALL_THRESHOLD)

    n_big = len(big_bars)
    cmap = plt.get_cmap("tab10" if n_big <= 10 else "tab20")

    bar_labels = [f"C{r}" for r, _ in big_bars] + ["Small\nclusters"]
    bar_values = [s for _, s in big_bars] + [small_total]
    bar_colors = [cmap(i / max(n_big - 1, 1)) for i in range(n_big)] + ["#b0b0b0"]

    fig, ax = plt.subplots(figsize=(max(7, len(bar_labels) * 0.8), 5))
    fig.patch.set_facecolor("#f9f9f9")
    ax.set_facecolor("#f9f9f9")

    bars = ax.bar(bar_labels, bar_values,
                  color=bar_colors, edgecolor="white", linewidth=0.8, width=0.65)

    for bar, val in zip(bars, bar_values):
        offset = bar.get_height() * 0.05
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + offset,
            f"{val:,}",
            ha="center", va="bottom", fontsize=8,
        )

    ax.set_ylabel("Number of conformers", fontsize=10)
    ax.set_title(
        f"ATP Conformer Cluster Sizes  ·  {len(clusters)} clusters total\n"
        f"'Small clusters' = {n_small_clusters} clusters with < {SMALL_THRESHOLD} members",
        fontsize=11,
    )
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{int(x):,}")
    )
    ax.set_ylim(bottom=1)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", default="output")
    parser.add_argument("--plot_dir",   default="pictures")
    parser.add_argument("--cutoff", type=float, default=None,
                        help="RMSD cutoff used during clustering (read from clusters.json if omitted)")
    args = parser.parse_args()
    d = args.output_dir
    Path(args.plot_dir).mkdir(parents=True, exist_ok=True)

    # ── Load ──────────────────────────────────────────────────────────────────
    print("Loading RMSD matrix …")
    rmsd_mat = np.load(f"{d}/rmsd_matrix.npy")
    with open(f"{d}/clusters.json") as f:
        meta = json.load(f)
    N, cutoff = rmsd_mat.shape[0], meta["cutoff_A"]

    # ── Re-cluster from saved matrix (fast) ───────────────────────────────────
    print("Re-clustering from saved distance matrix …")
    dist_list = [float(rmsd_mat[i, j]) for i in range(1, N) for j in range(i)]
    clusters = Butina.ClusterData(dist_list, N, cutoff, isDistData=True, reordering=True)
    clusters = sorted(clusters, key=len, reverse=True)
    print(f"  {len(clusters)} clusters recovered")

    labels = np.full(N, -1, dtype=int)
    for ci, cl in enumerate(clusters):
        for idx in cl:
            labels[idx] = ci

    big_indices = [i for i, cl in enumerate(clusters) if len(cl) >= SMALL_THRESHOLD]

    # ── Embedding ─────────────────────────────────────────────────────────────
    print("Embedding …")
    embedding = mds_embedding(rmsd_mat, f"{d}/embedding_mds.npy")

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("Plotting …")
    plot_embedding(embedding, labels, clusters, big_indices,
                   f"{args.plot_dir}/cluster_embedding_mds.png")
    plot_bar(clusters, f"{args.plot_dir}/cluster_sizes_bar.png")

    print(f"\nDone. Data in {d}/  |  Plots in {args.plot_dir}/")


if __name__ == "__main__":
    main()
