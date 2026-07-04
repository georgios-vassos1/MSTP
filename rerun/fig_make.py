#!/usr/bin/env python3
"""Regenerate paper figures Figure_3..Figure_5 from corrected-instance data (rerun/figdata/*.json).
Figure_1 (schematic) and Figure_2 (corrmat) are generated elsewhere.
Fonts ~10pt (slightly below the 12pt body text), vector PDF, clean layout.
Writes to ~/drayage/report_tlpr/figs/."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FD  = os.path.expanduser("~/drayage/MSTP/rerun/figdata")
OUT = os.path.expanduser("~/drayage/report_tlpr/figs")

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,          # body text is 12pt; figures slightly smaller but readable
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})
BLUE, GREY = "#2b6cb0", "#4a5568"

def load(name): return json.load(open(os.path.join(FD, name)))

# ---- Figure 3: sensitivity (lambda, rho, m), mean +/- sd ----
def fig2():
    lam = [(200,134379,33387),(700,380019,50766),(2000,1012886,103129)]
    rho = [(0.0,388819,63206),(0.2,370472,42446),(0.4,364373,41772)]
    mm  = [(1,400185,69989),(2,388355,65115),(4,388688,60783)]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    for a,(dat,xl,ttl) in zip(ax, [(lam,r"$\lambda$","(a) Demand rate"),
                                    (rho,r"$\rho_{\mathrm{cross}}$","(b) Correlation"),
                                    (mm,r"$m$","(c) Spot multiplier")]):
        xs=[d[0] for d in dat]; ys=[d[1]/1e3 for d in dat]; es=[d[2]/1e3 for d in dat]
        a.errorbar(range(len(xs)), ys, yerr=es, fmt="o-", color=BLUE, capsize=3, lw=1.5, ms=5)
        a.set_xticks(range(len(xs))); a.set_xticklabels(xs)
        a.set_xlabel(xl); a.set_title(ttl)
    ax[0].set_ylabel("Upper bound (×$10^3$)")
    fig.tight_layout(); fig.savefig(f"{OUT}/Figure_3.pdf", bbox_inches="tight"); plt.close(fig)

# (Removed from the paper: the regret-topo boxplot and the regret-vs-iters figure.
#  Regret across topologies is reported in the tables; SDDP convergence is now the
#  LB/UB bound-convergence figure below.)

# ---- Figure 5: capacity-opt convergence (obj, grad norm, step) ----
def fig5():
    d = load("capopt_traj.json"); it=d["iter"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax[0].plot(it, np.array(d["obj"])/1e3, "o-", color=BLUE, lw=1.5, ms=4); ax[0].set_title("(a) Objective"); ax[0].set_ylabel(r"reservation $+$ op. cost (×$10^3$)")
    ax[1].plot(it, d["gradnorm"], "o-", color=BLUE, lw=1.5, ms=4); ax[1].set_title(r"(b) $\|\nabla f\|_2$")
    ax[2].plot(it, d["step"], "o-", color=BLUE, lw=1.5, ms=4); ax[2].set_title(r"(c) $\|\Delta x\|_2$")
    for a in ax: a.set_xlabel("outer iteration")
    fig.tight_layout(); fig.savefig(f"{OUT}/Figure_4.pdf", bbox_inches="tight"); plt.close(fig)

# ---- Figure 4: SDDP convergence -- LB rises to a flat out-of-sample UB, gap closes ----
def fig6():
    d = load("convergence.json")
    it = np.array(d["iters"]); lb = np.array(d["lb"]); ub = np.array(d["ub"]); se = np.array(d["se"])
    m = it >= 5                                   # drop the iter=2 outlier for readability
    it, lb, ub, se = it[m], lb[m], ub[m], se[m]
    gap = 100 * (ub - lb) / lb
    GOLD = "#b7791f"
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    ax.fill_between(it, lb/1e3, ub/1e3, color=BLUE, alpha=0.12, lw=0, label="optimality gap")
    ax.fill_between(it, (ub-se)/1e3, (ub+se)/1e3, color=GREY, alpha=0.30, lw=0)
    ax.plot(it, ub/1e3, "s-", color=GREY, lw=1.3, ms=3.5, label="out-of-sample cost (UB)")
    ax.plot(it, lb/1e3, "o-", color=BLUE, lw=1.7, ms=4, label="lower bound (LB)")
    ax.set_xscale("log")
    ax.set_xticks([5, 10, 20, 50, 100, 300]); ax.set_xticklabels([5, 10, 20, 50, 100, 300])
    ax.set_xlabel("SDDP iterations"); ax.set_ylabel(r"policy cost (×$10^3$)")
    ax.legend(frameon=False, loc="center right", fontsize=8.5)
    ax2 = ax.twinx()
    ax2.plot(it, gap, ":", color=GOLD, lw=1.6, label="gap")
    ax2.set_ylabel("optimality gap (%)", color=GOLD); ax2.set_ylim(0, max(gap)*1.1)
    ax2.tick_params(axis="y", colors=GOLD); ax2.grid(False)
    ax2.spines["right"].set_visible(True); ax2.spines["right"].set_color(GOLD); ax2.spines["top"].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{OUT}/Figure_4.pdf", bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    fig2(); print("Figure_3 done")
    fig5(); print("Figure_5 done")
    fig6(); print("Figure_4 done")
    print("figures written to", OUT)
