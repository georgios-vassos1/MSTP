#!/usr/bin/env python3
"""Regenerate paper figures F2-F5 from the corrected-instance data (rerun/figdata/*.json).
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

# ---- F2: sensitivity (lambda, rho, m), mean +/- sd ----
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
    fig.tight_layout(); fig.savefig(f"{OUT}/F2_sensitivity.pdf", bbox_inches="tight"); plt.close(fig)

# ---- F3: regret distributions across topologies (near zero) ----
def fig3():
    d = load("regret_topo.json"); order=["t2x1","t1x2","t2x2"]; labels=[r"$2\times1$",r"$1\times2$",r"$2\times2$"]
    data=[np.array(d[k]) for k in order]
    fig, a = plt.subplots(figsize=(4.2, 2.8))
    bp=a.boxplot(data, tick_labels=labels, showfliers=True, widths=0.5,
                 patch_artist=True, medianprops=dict(color=GREY, lw=1.2),
                 flierprops=dict(marker='.', ms=3, mfc=GREY, mec=GREY, alpha=0.4))
    for patch in bp['boxes']: patch.set(facecolor=BLUE, alpha=0.35, lw=1.0)
    a.set_ylabel("Perfect-foresight regret (\\%)"); a.set_xlabel("Topology")
    a.axhline(0, color=GREY, lw=0.6, ls="--")
    fig.tight_layout(); fig.savefig(f"{OUT}/F3_regret_topology.pdf", bbox_inches="tight"); plt.close(fig)

# ---- F4: regret vs iteration budget (mean + IQR band) ----
def fig4():
    d = load("regret_iters.json"); its=sorted(int(k) for k in d)
    med=[np.median(d[str(i)]) for i in its]
    q1 =[np.percentile(d[str(i)],25) for i in its]
    q3 =[np.percentile(d[str(i)],75) for i in its]
    mean=[np.mean(d[str(i)]) for i in its]
    fig, a = plt.subplots(figsize=(4.4, 2.8))
    a.fill_between(its, q1, q3, color=BLUE, alpha=0.2, label="IQR")
    a.plot(its, med, "o-", color=BLUE, lw=1.5, ms=4, label="median")
    a.plot(its, mean, "s--", color=GREY, lw=1.2, ms=4, label="mean")
    a.set_xlabel("SDDP iterations"); a.set_ylabel("Perfect-foresight regret (\\%)")
    a.legend(frameon=False); a.set_ylim(bottom=0)
    fig.tight_layout(); fig.savefig(f"{OUT}/F4_regret_iterations.pdf", bbox_inches="tight"); plt.close(fig)

# ---- F5: capacity-opt convergence (obj, grad norm, step) ----
def fig5():
    d = load("capopt_traj.json"); it=d["iter"]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    ax[0].plot(it, np.array(d["obj"])/1e3, "o-", color=BLUE, lw=1.5, ms=4); ax[0].set_title("(a) Objective"); ax[0].set_ylabel(r"reservation $+$ op. cost (×$10^3$)")
    ax[1].plot(it, d["gradnorm"], "o-", color=BLUE, lw=1.5, ms=4); ax[1].set_title(r"(b) $\|\nabla f\|_2$")
    ax[2].plot(it, d["step"], "o-", color=BLUE, lw=1.5, ms=4); ax[2].set_title(r"(c) $\|\Delta x\|_2$")
    for a in ax: a.set_xlabel("outer iteration")
    fig.tight_layout(); fig.savefig(f"{OUT}/F5_cap_opt_convergence.pdf", bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    fig2(); print("F2 done")
    fig3(); print("F3 done")
    fig4(); print("F4 done")
    fig5(); print("F5 done")
    print("figures written to", OUT)
