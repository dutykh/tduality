#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the pair energy and the running interaction dimension on the square torus.

The data are the reference outputs stored in convergence/, which the two Julia scripts of codes/
write again. From this folder,
    python3 plot_readme_figure.py
writes tduality-torus.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "convergence"
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

RATIOS = ((0.02, "0p020"), (0.05, "0p050"), (0.10, "0p100"))

hs.configure(family="latin-modern")
plt.rcParams["savefig.dpi"] = 265
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.7), constrained_layout=True)
handles, labels = [], []
for i, (a, tag) in enumerate(RATIOS):
    style = hs.series(i)
    colour, marker = style["color"], style["marker"]
    xi, torus, plane = np.loadtxt(DATA / f"torus_potential_a_{tag}.csv", delimiter=",",
                                  skiprows=1, unpack=True)
    left.plot(xi, torus, color=colour, lw=1.1, marker=marker, markevery=(25, 50), markersize=3.6,
              markerfacecolor="white", markeredgecolor=colour)
    left.plot(xi, plane, color=colour, lw=0.9, linestyle=hs.DASHES[1])
    xi, dim_torus, dim_plane = np.loadtxt(DATA / f"running_dimension_a_{tag}.csv", delimiter=",",
                                          skiprows=1, unpack=True)
    right.plot(xi, dim_torus, color=colour, lw=1.1, marker=marker, markevery=(25, 50), markersize=3.6,
               markerfacecolor="white", markeredgecolor=colour)
    right.plot(xi, dim_plane, color=colour, lw=0.9, linestyle=hs.DASHES[1])
    handles.append(Line2D([], [], color=colour, lw=1.1, marker=marker, markersize=3.6,
                          markerfacecolor="white", markeredgecolor=colour))
    labels.append(rf"$\ell_0/L = {a:.2f}$")

handles += [Line2D([], [], color=hs.INK, lw=1.1), Line2D([], [], color=hs.INK, lw=0.9,
                                                         linestyle=hs.DASHES[1])]
labels += ["square torus", "infinite plane"]

for axis in (left, right):
    axis.set_xlim(0.0, 0.5)
    axis.set_xlabel(r"separation $\xi = r/L$")
left.set_ylim(0.0, 0.55)
left.set_ylabel(r"pair energy $\Delta V/Q^2$")
hs.panel(left, "(a) pair energy along the axis")
right.set_ylim(0.9, 3.1)
right.set_yticks([1.0, 1.5, 2.0, 2.5, 3.0])
right.set_ylabel(r"running dimension $D_\Delta$")
hs.panel(right, "(b) running interaction dimension")

hs.legend_outside(fig, where="lower", ncol=5, handles=handles, labels=labels)
hs.save(fig, HERE / "tduality-torus.png")
print("written", HERE / "tduality-torus.png")
