import numpy as np
import h5py
from matplotlib import pyplot as plt
import os

# ***** Parameters *****
filename = "../Data.h5"
outdir   = "Frames"

ymin = -0.5
ymax = 14. 

figsize = (10., 7.)
dpi     = 200

time_x =  0.65
time_y = -0.08
# **********************


##if os.path.isdir(outdir):
##    raise RuntimeError(f"Directory '{outdir}' already exists. Please remove it or change the 'outdir' variable before re-running this script.")
os.makedirs(outdir, exist_ok = False)

with h5py.File(filename, "r") as f:
    t         = f["Time"][()]
    x         = f["Position"][()]
    out_every = f["Output frequency"][()]

    NT = len(t)
    ##NX = len(x)

    assert len(out_every) == 1
    out_every = out_every[0]

    for n in range(NT):
        it = n*out_every
        u  = f[f"Function/Iteration {it}"][()]

        fig = plt.figure(figsize = figsize, dpi = dpi)
        ax  = fig.gca()
        plt.xlabel("$x$", fontsize = 12.)
        plt.ylabel(r"$u\left(x\right)$", fontsize = 12.)
        plt.grid(ls = "--", lw = 0.5, alpha = 0.5)
        plt.text(time_x, time_y, f"t = {t[n]:.2e}", transform = ax.transAxes, fontsize = 12., fontweight = "bold")
        plt.ylim(ymin, ymax)
        plt.plot(x, u, ls = "-",  lw = 1., marker = ".", markersize = 2., color = "dodgerblue")
        plt.tight_layout()
        plt.savefig(f"{outdir}/it_{n:06d}.png")
        plt.close()

        print(f"Iteration {it}/{NT*out_every}, time {t[n]:.2e} done")
