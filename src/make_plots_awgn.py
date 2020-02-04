"""
A Python script to make beautiful plots
"""

#%% Imports
from contextlib import suppress
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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
results_awgn = loadmat("../results/resultsawgn.mat")

tebnos = results_awgn["tebnos"][0]

awgn_pskMs = results_awgn["pskMs"][0]
awgn_pskMs.sort()
awgn_respsk_teo = results_awgn["resultspsk"][1, ...]
awgn_respsk = results_awgn["resultspsk"][0, ...]

awgn_qamMs = results_awgn["qamMs"][0]
awgn_qamMs.sort()
awgn_resqam_teo = results_awgn["resultsqam"][1, ...]
awgn_resqam = results_awgn["resultsqam"][0, ...]

awgn_fskMs = results_awgn["fskMs"][0]
awgn_fskMs.sort()
awgn_resfsk_teo = results_awgn["resultsfsk"][1, ...]
awgn_resfsk = results_awgn["resultsfsk"][0, ...]

awgn_respsk_cyclic = results_awgn["resultspsk_crc"]
awgn_resqam_cyclic = results_awgn["resultsqam_crc"]
awgn_resfsk_cyclic = results_awgn["resultsfsk_crc"]

awgn_cyclicCodes = [
    list(results_awgn["cyclicCodes"][0][i][0])
    for i in range(len(results_awgn["cyclicCodes"][0]))
]

# Markers to use on the cyclic FEC plots
awgn_cyclic_markers = ["^", "v"]

awgn_respsk_conv = results_awgn["resultspsk_conv"]
awgn_resqam_conv = results_awgn["resultsqam_conv"]
awgn_resfsk_conv = results_awgn["resultsfsk_conv"]
awgn_trellis_args = ["2/3 Feedforward", "1/2 Feedforward", "1/2 Feedback"]

awgn_trellis_markers = ["o", "D", "s"]


#%% AWGN PSK results
plt.figure(1)
plt.title("BER vs $E_b/N_o$ for PSK modulation")
artists = list()
for (idx, M) in enumerate(awgn_pskMs):
    plt.plot(tebnos, awgn_respsk_teo[idx], label=f"M = {M}", color=colors[idx])
    plt.plot(
        tebnos,
        awgn_respsk[idx],
        label=f"M = {M}",
        color=colors[idx],
        marker="*",
        linestyle="None",
    )
    artists.append(
        plt.Line2D((0, 1), (0, 0), color=colors[idx], marker="*", linestyle="-")
    )
plt.yscale("log", nonposy="mask")
ax = plt.gca()
ax.set_ylim(bottom=10e-5, top=1)
ax.set_xlim(left=min(tebnos), right=max(tebnos))
ax.grid(True, which="both")
# Get artists and labels for legend and chose which ones to display
handles, labels = ax.get_legend_handles_labels()
# Create custom artists
simArtist = plt.Line2D((0, 1), (0, 0), color="k", marker="*", linestyle="")
anyArtist = plt.Line2D((0, 1), (0, 0), color="k")
# Create legend from custom artist/label lists
todisplay = range(0, len(awgn_pskMs) * 2, 2)
ax.legend(
    artists + [simArtist, anyArtist],
    [label for i, label in enumerate(labels) if i in todisplay]
    + ["Simulation", "Analytic"],
    framealpha=1,
)
ax.set_xlabel("$E_b/N_o$ (dB)")
ax.set_ylabel("BER")
ax.xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))

plt.savefig("../results/figure_1.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_1.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_1.pdf",  bbox_inches="tight")


#%% AWGN QAM results
plt.figure(2)
plt.title("BER vs $E_b/N_o$ for QAM modulation")
artists = list()
for (idx, M) in enumerate(awgn_qamMs):
    plt.plot(tebnos, awgn_resqam_teo[idx], label=f"M = {M}", color=colors[idx])
    plt.plot(
        tebnos,
        awgn_resqam[idx],
        label=f"M = {M}",
        color=colors[idx],
        marker="*",
        linestyle="None",
    )
    artists.append(
        plt.Line2D((0, 1), (0, 0), color=colors[idx], marker="*", linestyle="-")
    )
ax = plt.gca()
# Get artists and labels for legend and chose which ones to display
handles, labels = ax.get_legend_handles_labels()
# Create custom artists
simArtist = plt.Line2D((0, 1), (0, 0), color="k", marker="*", linestyle="")
anyArtist = plt.Line2D((0, 1), (0, 0), color="k")
# Create legend from custom artist/label lists
todisplay = range(0, len(awgn_qamMs) * 2, 2)
ax.legend(
    artists + [simArtist, anyArtist],
    [label for i, label in enumerate(labels) if i in todisplay]
    + ["Simulation", "Analytic"],
    framealpha=1,
)
plt.yscale("log", nonposy="mask")
ax.set_ylim(bottom=10e-5, top=1)
ax.set_xlim(left=min(tebnos), right=max(tebnos))
ax.grid(True, which="both")
ax.set_xlabel("$E_b/N_o$ (dB)")
ax.set_ylabel("BER")
ax.xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))

plt.savefig("../results/figure_2.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_2.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_2.pdf",  bbox_inches="tight")

#%% AWGN FSK results
plt.figure(3)
plt.title("BER vs $E_b/N_o$ for FSK modulation")
artists = list()
for (idx, M) in enumerate(awgn_fskMs):
    plt.plot(tebnos, awgn_resfsk_teo[idx], label=f"M = {M}", color=colors[idx])
    plt.plot(
        tebnos,
        awgn_resfsk[idx],
        label=f"M = {M}",
        color=colors[idx],
        marker="*",
        linestyle="None",
    )
    artists.append(
        plt.Line2D((0, 1), (0, 0), color=colors[idx], marker="*", linestyle="-")
    )
ax = plt.gca()
# Get artists and labels for legend and chose which ones to display
handles, labels = ax.get_legend_handles_labels()
# Create custom artists
simArtist = plt.Line2D((0, 1), (0, 0), color="k", marker="*", linestyle="")
anyArtist = plt.Line2D((0, 1), (0, 0), color="k")
# Create legend from custom artist/label lists
todisplay = range(0, len(awgn_fskMs) * 2, 2)
ax.legend(
    artists + [simArtist, anyArtist],
    [label for i, label in enumerate(labels) if i in todisplay]
    + ["Simulation", "Analytic"],
    framealpha=1,
)
plt.yscale("log", nonposy="mask")
ax.set_ylim(bottom=10e-5, top=1)
ax.set_xlim(left=min(tebnos), right=max(tebnos))
ax.grid(True, which="both")
ax.set_xlabel("$E_b/N_o$ (dB)")
ax.set_ylabel("BER")
ax.xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))

plt.savefig("../results/figure_3.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_3.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_3.pdf",  bbox_inches="tight")

#%% AWG PSK Cyclic FEC results
fig, axs = plt.subplots(
    len(awgn_pskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=4,
    figsize=(6, 10),
)

for (idx, M) in enumerate(awgn_pskMs):
    axs[idx].plot(tebnos, awgn_respsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_respsk_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"n = {n}, k = {k}",
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left")
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for PSK modulation w/ cyclic FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_4.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_4.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_4.pdf",  bbox_inches="tight")

#%% AWG QAM Cyclic FEC results
fig, axs = plt.subplots(
    len(awgn_qamMs), sharex=True, gridspec_kw={"wspace": 0, "hspace": 0.1}, num=5
)

for (idx, M) in enumerate(awgn_qamMs):
    axs[idx].plot(tebnos, awgn_resqam_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_resqam_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"n = {n}, k = {k}",
        )
    axs[idx].legend(framealpha=1)
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for QAM modulation w/ cyclic FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_5.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_5.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_5.pdf",  bbox_inches="tight")


#%% AWG FSK Cyclic FEC results
fig, axs = plt.subplots(
    len(awgn_fskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=6,
    figsize=(6, 10),
)

for (idx, M) in enumerate(awgn_fskMs):
    axs[idx].plot(tebnos, awgn_resfsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_resfsk_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"n = {n}, k = {k}",
        )
    axs[idx].legend(framealpha=1, loc="lower left")
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for FSK modulation w/ cyclic FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_6.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_6.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_6.pdf",  bbox_inches="tight")


#%% AWG PSK Convolutional FEC results
fig, axs = plt.subplots(
    len(awgn_pskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=7,
    figsize=(6, 10),
)

for (idx, M) in enumerate(awgn_pskMs):
    axs[idx].plot(tebnos, awgn_respsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_respsk_conv[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label=awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left")
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-8, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for PSK modulation w/ convolutional FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_7.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_7.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_7.pdf",  bbox_inches="tight")


#%% AWG QAM Convolutional FEC results
fig, axs = plt.subplots(
    len(awgn_qamMs), sharex=True, gridspec_kw={"wspace": 0, "hspace": 0.1}, num=8
)

for (idx, M) in enumerate(awgn_qamMs):
    axs[idx].plot(tebnos, awgn_resqam_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_resqam_conv[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label=awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left")
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-8, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for QAM modulation w/ convolutional FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_8.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_8.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_8.pdf",  bbox_inches="tight")


#%% AWG FSK Convolutional FEC results
fig, axs = plt.subplots(
    len(awgn_fskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=9,
    figsize=(6, 10),
)

for (idx, M) in enumerate(awgn_fskMs):
    axs[idx].plot(tebnos, awgn_resfsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_resfsk_conv[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label=awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left")
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-8, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for FSK modulation w/ convolutional FEC")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_9.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_9.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_9.pdf",  bbox_inches="tight")


#%% Plot common elements
awgn_commonMs = list(
    set().union(
        set(awgn_pskMs).intersection(awgn_qamMs),
        set(awgn_pskMs).intersection(awgn_fskMs),
        set(awgn_fskMs).intersection(awgn_qamMs),
    )
)
awgn_commonMs.sort()
fig, axs = plt.subplots(
    len(awgn_commonMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=10,
    figsize=(6, 10),
)

for idx, M in enumerate(awgn_commonMs):
    simArtist = plt.Line2D((0, 1), (0, 0), color="k", marker="*", linestyle="")
    anyArtist = plt.Line2D((0, 1), (0, 0), color="k")
    artists = list()
    labels = list()
    try:
        idxpsk = list(awgn_pskMs).index(M)
        axs[idx].plot(
            tebnos, awgn_respsk_teo[idxpsk], label=f"PSK, M = {M}", color=colors[0]
        )
        axs[idx].plot(
            tebnos,
            awgn_respsk[idxpsk],
            label=f"M = {M}",
            color=colors[0],
            marker="*",
            linestyle="None",
        )
        artists.append(
            plt.Line2D((0, 1), (0, 0), color=colors[0], marker="*", linestyle="-")
        )
        labels.append(f"PSK, M = {M}")
    except ValueError:
        pass
    try:
        idxqam = list(awgn_qamMs).index(M)
        axs[idx].plot(
            tebnos, awgn_resqam_teo[idxqam], label=f"QAM, M = {M}", color=colors[2]
        )
        axs[idx].plot(
            tebnos,
            awgn_resqam[idxqam],
            label=f"M = {M}",
            color=colors[2],
            marker="*",
            linestyle="None",
        )
        artists.append(
            plt.Line2D((0, 1), (0, 0), color=colors[2], marker="*", linestyle="-")
        )
        labels.append(f"QAM, M = {M}")
    except ValueError:
        pass
    try:
        idxfsk = list(awgn_fskMs).index(M)
        axs[idx].plot(
            tebnos, awgn_resfsk_teo[idxfsk], label=f"FSK, M = {M}", color=colors[3]
        )
        axs[idx].plot(
            tebnos,
            awgn_resfsk[idxfsk],
            label=f"M = {M}",
            color=colors[3],
            marker="*",
            linestyle="None",
        )
        artists.append(
            plt.Line2D((0, 1), (0, 0), color=colors[3], marker="*", linestyle="-")
        )
        labels.append(f"FSK, M = {M}")
    except ValueError:
        pass
    axs[idx].legend(
        artists + [simArtist, anyArtist],
        labels + ["Simulation", "Analytic"],
        framealpha=1,
        fontsize=6,
        loc="lower left",
    )
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("Comparison of BER vs $E_b/N_o$ for different modulations")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_10.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_10.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_10.pdf",  bbox_inches="tight")


#%% Compare the FECs for PSK
fig, axs = plt.subplots(
    len(awgn_pskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=11,
    figsize=(6, 10),
)

for (idx, M) in enumerate(awgn_pskMs):
    axs[idx].plot(tebnos, awgn_respsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_respsk_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"Cyclic; n = {n}, k = {k}",
        )
    didx = len(awgn_cyclicCodes)
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_respsk_conv[cidx, idx],
            color=colors[didx + cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label="Conv.; " + awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left", ncol=2)
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for PSK mod. w/ cyclic or conv. FECs")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_11.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_11.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_11.pdf",  bbox_inches="tight")


#%% Compare the FECs for QAM
fig, axs = plt.subplots(
    len(awgn_qamMs), sharex=True, gridspec_kw={"wspace": 0, "hspace": 0.1}, num=12
)

for (idx, M) in enumerate(awgn_qamMs):
    axs[idx].plot(tebnos, awgn_resqam_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_resqam_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"Cyclic; n = {n}, k = {k}",
        )
    didx = len(awgn_cyclicCodes)
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_resqam_conv[cidx, idx],
            color=colors[didx + cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label="Conv.; " + awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left", ncol=2)
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for QAM mod. w/ cyclic or conv. FECs")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_12.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_12.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_12.pdf",  bbox_inches="tight")


#%% Compare the FECs for FSK
fig, axs = plt.subplots(
    len(awgn_fskMs),
    sharex=True,
    gridspec_kw={"wspace": 0, "hspace": 0.1},
    num=13,
    figsize=(6, 8),
)

for (idx, M) in enumerate(awgn_fskMs):
    axs[idx].plot(tebnos, awgn_resfsk_teo[idx], label=f"M = {M}", color=colors[0])
    for (cidx, (n, k)) in enumerate(awgn_cyclicCodes):
        axs[idx].plot(
            tebnos,
            awgn_resfsk_cyclic[cidx, idx],
            color=colors[cidx + 1],
            marker=awgn_cyclic_markers[cidx],
            linestyle="--",
            markersize=4,
            label=f"Cyclic; n = {n}, k = {k}",
        )
    didx = len(awgn_cyclicCodes)
    for (cidx, trelname) in enumerate(awgn_trellis_args):
        axs[idx].plot(
            tebnos,
            awgn_resfsk_conv[cidx, idx],
            color=colors[didx + cidx + 1],
            marker=awgn_trellis_markers[cidx],
            linestyle="--",
            markersize=4,
            label="Conv.; " + awgn_trellis_args[cidx],
        )
    axs[idx].legend(framealpha=1, fontsize=6, loc="lower left", ncol=2)
    axs[idx].set_yscale("log", nonposy="mask")
    axs[idx].set_ylim(bottom=10e-6, top=1)
    axs[idx].grid(True, which="both")
    axs[idx].xaxis.set_major_locator(ticker.MultipleLocator(awgn_xtickspacing))
axs[0].set_title("BER vs $E_b/N_o$ for FSK mod. w/ cyclic or conv. FECs")
axs[-1].set_xlabel("$E_b/N_o$ (dB)")
axs[-1].set_xlim(left=min(tebnos), right=max(tebnos))

plt.savefig("../results/figure_13.png", dpi=600, bbox_inches="tight")
plt.savefig("../report/figures/pgf/figure_13.pgf",  bbox_inches="tight")
plt.savefig("../report/figures/pdf/figure_13.pdf",  bbox_inches="tight")


#%% Show plots
plt.show()

#%%
