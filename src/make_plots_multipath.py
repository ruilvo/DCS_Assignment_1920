"""
A Python script to make beautiful plots
"""

#%% Imports
from contextlib import suppress
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import matplotlib


#%% Import ipython and adjust for qt5 backend, keeping compatibility with
# non-ipython envs
with suppress(ImportError):
    from IPython import get_ipython
with suppress(AttributeError, NameError):
    # List available APIs
    get_ipython().run_line_magic("matplotlib", "-l")
    get_ipython().run_line_magic("matplotlib", "qt5")

print(plt.get_backend())
#%% Matplotlib adjustments
# Adjust matplotlib qt5 backend for high DPI monitors
if plt.get_backend() == "Qt5Agg":
    from matplotlib.backends.qt_compat import QtWidgets  # pylint: disable=C0412
    from sys import argv

    qApp = QtWidgets.QApplication(argv)
    plt.matplotlib.rcParams["figure.dpi"] = qApp.desktop().physicalDpiX()

# Set LaTeX compatibility settings
# matplotlib.use("pgf")
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "sans-serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
    }
)

plt.rcParams["figure.autolayout"] = True

# Adjust figure default size
plt.rcParams["figure.figsize"] = (6, 4)  # Default fig size (w, h))

# Get the default color scheme
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

awgn_xtickspacing = 5


#%% Get AWGN results from files
# Load the file
results_mpath = loadmat("../results/resultsmultipath.mat")

pskMs = results_mpath["pskMs"][0]
pskMs.sort()

qamMs = results_mpath["qamMs"][0]
qamMs.sort()

tdels = results_mpath["timedelays"][0] / 1e-6
rpwrs = results_mpath["relpowers"][0]
npaths = results_mpath["npaths"][0]

res_psk_delays = results_mpath["results_psk_delays"]
res_psk_npaths = results_mpath["results_psk_npaths"]
res_psk_pwrs = results_mpath["results_psk_powers"]

res_qam_delays = results_mpath["results_qam_delays"]
res_qam_npaths = results_mpath["results_qam_npaths"]
res_qam_pwrs = results_mpath["results_qam_powers"]


#%% Load multipath results

results_ofdm = loadmat("../results/resultsmultipathofdm.mat")

pskMs_ofdm = results_ofdm["pskMs"][0]
pskMs_ofdm.sort()

qamMs_ofdm = results_ofdm["qamMs"][0]
qamMs_ofdm.sort()

tdels_ofdm = results_ofdm["timedelays"][0] / 1e-6
rpwrs_ofdm = results_ofdm["relpowers"][0]
npaths_ofdm = results_ofdm["npaths"][0]

res_psk_delays_ofdm = results_ofdm["results_psk_delays"]
res_psk_npaths_ofdm = results_ofdm["results_psk_npaths"]
res_psk_pwrs_ofdm = results_ofdm["results_psk_powers"]

res_qam_delays_ofdm = results_ofdm["results_qam_delays"]
res_qam_npaths_ofdm = results_ofdm["results_qam_npaths"]
res_qam_pwrs_ofdm = results_ofdm["results_qam_powers"]


#%% Multipath on normal modulation
plt.figure(1)
plt.title("BER vs excess channel delay")
for idx, M in enumerate(pskMs):
    plt.plot(
        tdels,
        res_psk_delays[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
didx = len(pskMs)
for idx, M in enumerate(qamMs):
    plt.plot(
        tdels,
        res_qam_delays[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
plt.yscale("log", nonposy="clip")
ax = plt.gca()
ax.set_xlabel("Excess delay ($\\mu$s)")
ax.set_ylim(top=0.1, bottom=10e-6)
ax.set_xlim(left=min(tdels), right=max(tdels))
ax.legend(framealpha=1, fontsize=6, loc="lower right")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_1.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_1.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_1.pdf", bbox_inches="tight")


#%% Multipath on OFDM modulation
plt.figure(2)
plt.title("BER vs excess channel delay for OFDM modulation")
for idx, M in enumerate(pskMs_ofdm):
    plt.plot(
        tdels_ofdm,
        res_psk_delays_ofdm[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
didx = len(pskMs_ofdm)
for idx, M in enumerate(qamMs_ofdm):
    plt.plot(
        tdels_ofdm,
        res_qam_delays_ofdm[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
plt.yscale("log", nonposy="clip")
ax = plt.gca()
ax.set_xlabel("Excess delay ($\\mu$s)")
ax.set_xlim(left=min(tdels_ofdm), right=max(tdels_ofdm))
ax.set_ylim(top=0.1, bottom=10e-6)
ax.legend(framealpha=1, fontsize=6, loc="upper left")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_2.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_2.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_2.pdf", bbox_inches="tight")


#%% Compare for PSK and QAM
common_pskMs = list(set(pskMs).intersection(pskMs_ofdm))
common_qamMs = list(set(qamMs).intersection(qamMs_ofdm))
fig = plt.figure(constrained_layout=False, num=3)
fig.suptitle(
    'BER for PSK (left) and QAM (right) for both "classical" and\n'
    + "OFDM modulation scenarios",
    y=0.98,
)
nrows = len(common_qamMs) * len(common_pskMs) + 1
gs = GridSpec(
    nrows=nrows,
    ncols=2,
    figure=fig,
    top=0.88,
    left=0.05,
    right=0.95,
    bottom=0.08,
    wspace=0.05,
)
# Now we plot the stuff
axs = list()
for idx, M in enumerate(common_pskMs):
    # This seems counter intuitive, but we need to plot in steps of the OTHER
    ax = fig.add_subplot(
        gs[idx * len(common_qamMs) : idx * len(common_qamMs) + len(common_qamMs), 0]
    )
    axs.append(ax)
    ax.plot(
        tdels,
        res_psk_delays[idx],
        label=f"without OFDM",
        color=colors[2],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
    ax.plot(
        tdels_ofdm,
        res_psk_delays_ofdm[idx],
        label=f"with OFDM",
        color=colors[0],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
    plt.yscale("log", nonposy="clip")
    ax.set_xlim(left=0, right=max(tdels_ofdm))
    toplim = max((res_psk_delays[idx].max(), res_psk_delays_ofdm[idx].max()))
    ax.set_ylim(top=toplim, bottom=10e-6)
    # Turn off ticks, expect for the bottom one
    if idx != (len(common_pskMs) - 1):
        ax.set_xticklabels([])
        ax.set_xticks([])
    else:
        ax.set_xlabel("Excess delay ($\\mu$s)")
        # make the x ticks integers, not floats
        xint = []
        locs = ax.get_xticks()
        for each in locs:
            xint.append(int(each))
            ax.set_xticks(xint)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.text(
        0.1,
        toplim * 0.05,
        f"PSK; M = {M}",
        bbox=dict(boxstyle="square", facecolor="white"),
        fontsize=6,
    )
    ax.axvline(0.8, color="red")
    ax.legend(framealpha=1, fontsize=6, loc="lower right", ncol=1)

for idx, M in enumerate(common_qamMs):
    # This seems counter intuitive, but we need to plot in steps of the OTHER
    ax = fig.add_subplot(
        gs[idx * len(common_pskMs) : idx * len(common_pskMs) + len(common_pskMs), 1]
    )
    ax.plot(
        tdels,
        res_qam_delays[0],
        label=f"without OFDM",
        color=colors[2],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
    ax.plot(
        tdels_ofdm,
        res_qam_delays_ofdm[idx],
        label=f"with OFDM",
        color=colors[0],
        marker="None",
        linestyle="-",
        linewidth=1,
    )
    plt.yscale("log", nonposy="clip")
    ax.set_xlim(left=0, right=max(tdels_ofdm))
    toplim = max((res_qam_delays[idx].max(), res_qam_delays_ofdm[idx].max())) + 0.003
    ax.set_ylim(top=toplim, bottom=10e-6)
    ax.text(
        0.1,
        toplim * 0.05,
        f"QAM; M = {M}",
        bbox=dict(boxstyle="square", facecolor="white"),
        fontsize=6,
    )
    # Turn off ticks, expect for the bottom one
    if idx != (len(common_qamMs) - 1):
        ax.set_xticklabels([])
        ax.set_xticks([])
    else:
        ax.set_xlabel("Excess delay ($\\mu$s)")
        # make the x ticks integers, not floats
        xint = []
        locs = ax.get_xticks()
        for each in locs:
            xint.append(int(each))
            ax.set_xticks(xint)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.axvline(0.8, color="red")
    ax.legend(framealpha=1, fontsize=6, loc="lower right", ncol=1)

plt.savefig("../results/figure_mp_3.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_3.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_3.pdf", bbox_inches="tight")


#%% Performance regarding number of paths
plt.figure(4)
plt.title("BER vs number of paths")
for idx, M in enumerate(pskMs):
    plt.plot(
        npaths,
        res_psk_npaths[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="^",
        linestyle="--",
        # linewidth=1,
    )
didx = len(pskMs)
for idx, M in enumerate(qamMs):
    plt.plot(
        npaths,
        res_qam_npaths[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="v",
        linestyle="--",
        # linewidth=1,
    )
plt.yscale("log", nonposy="mask")
ax = plt.gca()
ax.set_xlabel("Number of different paths")
ax.set_ylim(top=0.5, bottom=10e-5)
ax.set_xlim(left=0, right=max(npaths))
start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 2))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.legend(framealpha=1, fontsize=8, loc="lower right")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_4.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_4.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_4.pdf", bbox_inches="tight")


#%% Performance regarding rel. chan. power
plt.figure(5)
plt.title("BER vs relative chan. power")
for idx, M in enumerate(pskMs):
    plt.plot(
        rpwrs,
        res_psk_pwrs[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="^",
        linestyle="--",
        # linewidth=1,
    )
didx = len(pskMs)
for idx, M in enumerate(qamMs):
    plt.plot(
        rpwrs,
        res_qam_pwrs[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="v",
        linestyle="--",
        # linewidth=1,
    )
plt.yscale("log", nonposy="mask")
ax = plt.gca()
ax.set_xlabel("Second chan. power (dB to first chan)")
ax.set_ylim(top=0.5, bottom=10e-6)
ax.set_xlim(left=-15, right=5)
ax.legend(framealpha=1, fontsize=8, loc="lower right")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_5.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_5.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_5.pdf", bbox_inches="tight")


#%% Performance regarding number of paths
plt.figure(6)
plt.title("BER vs number of paths for OFDM")
for idx, M in enumerate(pskMs):
    plt.plot(
        npaths_ofdm,
        res_psk_npaths_ofdm[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="^",
        linestyle="--",
        # linewidth=1,
    )
didx = len(pskMs)
for idx, M in enumerate(qamMs):
    plt.plot(
        npaths_ofdm,
        res_qam_npaths_ofdm[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="v",
        linestyle="--",
        # linewidth=1,
    )
plt.yscale("log", nonposy="mask")
ax = plt.gca()
ax.set_xlabel("Number of different paths")
ax.set_ylim(top=0.5, bottom=10e-5)
ax.set_xlim(left=0, right=max(npaths))
start, end = ax.get_xlim()
ax.xaxis.set_ticks(np.arange(start, end, 2))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.legend(framealpha=1, fontsize=8, loc="lower right")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_6.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_6.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_6.pdf", bbox_inches="tight")


#%% Performance regarding rel. chan. power for OFDM
plt.figure(7)
plt.title("BER vs relative chan. power for OFDM")
for idx, M in enumerate(pskMs):
    plt.plot(
        rpwrs_ofdm,
        res_psk_pwrs_ofdm[idx],
        label=f"PSK; M = {M}",
        color=colors[idx],
        marker="^",
        linestyle="--",
        # linewidth=1,
    )
didx = len(pskMs)
for idx, M in enumerate(qamMs):
    plt.plot(
        rpwrs_ofdm,
        res_qam_pwrs_ofdm[idx],
        label=f"QAM; M = {M}",
        color=colors[didx + idx],
        marker="v",
        linestyle="--",
        # linewidth=1,
    )
plt.yscale("log", nonposy="mask")
ax = plt.gca()
ax.set_xlabel("Second chan. power (dB to first chan)")
ax.set_ylim(top=0.5, bottom=10e-6)
ax.set_xlim(left=-15, right=5)
ax.legend(framealpha=1, fontsize=8, loc="lower right")
ax.grid(True, which="both")
ax.set_ylabel("BER")

plt.savefig("../results/figure_mp_7.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_mp_7.pgf", bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_mp_7.pdf", bbox_inches="tight")

#%%
plt.show()
