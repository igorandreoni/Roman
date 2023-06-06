import glob
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib.ticker import MaxNLocator

lines_dict = {"Rubin": "--", "Roman":"-"}

title = r"GW170817 best fit at 1Gpc, $\theta$=26 deg"
out_filename = "figures_lc/prelim_lc_polar.pdf"
models = ["/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.005-0.20-0.20_wind0.050-0.05_theta25.84_z0.2.dat"]
#models = ["/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.005-0.20-0.20_wind0.050-0.05_theta78.46_z0.2.dat"]

if __name__ == "__main__":
    allfiles = models
    for filename in allfiles:
        # Read the data file
        t = ascii.read(filename)
        fig, ax = plt.subplots()
        # Create a zoomed inset.
        axins = ax.inset_axes([0.51, 0.51, 0.48, 0.48])
        n = len(t.colnames) - 1
        ax.set_prop_cycle('color',[plt.cm.jet(i) for i in np.linspace(0, 1, n)])
        for filt in t.colnames:
            # Ignore time column
            if filt == "t[days]":
                continue
            if "lsst" in filt:
                linestyle = lines_dict["Rubin"]
            else:
                linestyle = lines_dict["Roman"]
            # Interpolate to reduce noise
            idx = [i for i in np.arange(len(t)) if np.isnan(t[filt][i]) == False]
            idx = [i for i in idx if not t[filt][i] == np.inf]
            xnew = np.linspace(t["t[days]"][idx].min(), t["t[days]"][idx].max(), len(idx))
            f = interpolate.interp1d(t["t[days]"][idx],  t[filt][idx])
            if filt[0] == "l":
                label = f"{filt.replace('lsst', 'LSST ')}"
            else:
                label = f"Roman {filt.replace('f', 'F')}"
            ax.plot(xnew, f(xnew), linestyle=linestyle, label=label)
            # Plot also in the inset
            axins.plot(xnew, f(xnew), linestyle=linestyle)
        # Colors
        #from matplotlib.cm import get_cmap
        #name = "Oranges"
        #ax.cmap=mpl.colormaps[name]

        # Plot Roman sensitivity
        ax.plot([0, 55], [24.03,24.03], color='k', linewidth=3, linestyle='--', label="LSST r limit 30s")
        ax.plot([0, 55], [25.0,25.0], color='orchid', linewidth=3, linestyle='-', label="F129 limit 55s")
        ax.plot([0, 55], [28.0,28.0], color='fuchsia', linewidth=3, linestyle='-', label="F129 limit 1hr")
        # And in the insets
        axins.plot([0, 55], [24.03,24.03], color='k', linewidth=3, linestyle='--', label="LSST r limit 30s")
        axins.plot([0, 55], [25.0,25.0], color='orchid', linewidth=3, linestyle='-', label="F129 limit 55s")

        # Set plot parameters
        plt.rcParams['xtick.labelsize']=11
        ax.legend(fontsize=11, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_xlim([0.5, 21])
        ax.set_ylim([29, 22.5])
        ax.set_xlabel("Time from merger (days)", fontsize=13)
        ax.set_ylabel("AB Magnitude", fontsize=13)
        ax.set_title(title, fontsize=13)
        ax.tick_params(axis='both', labelsize=13)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        # subregion of the original image
        x1, x2, y1, y2 = 0.5, 3.5, 26, 22.5
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        #axins.set_xticklabels([])
        axins.set_yticklabels([])
        #ax.indicate_inset_zoom(axins, edgecolor="black")

        plt.savefig(out_filename, bbox_inches='tight')
        plt.show()

