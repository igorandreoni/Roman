import glob

import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib.ticker import MaxNLocator
from astropy.cosmology import Planck18 as cosmo
import astropy.cosmology as ac
import astropy.units as u

########

def chart(filters, days, results, exp, z):
    import math
    from collections import OrderedDict
    from matplotlib import colormaps
    from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator, ScalarFormatter)

    # for the filter colors:
    cmap = colormaps['jet']
    filt_color = {}
    for i, ff in zip(np.linspace(0, 1, len(filters)), filters):
        filt_color[ff] = cmap(i)
    ylabels = list(results.keys())
    ilen=len(ylabels)
    pos = np.arange(0.5,ilen*0.5+0.5,0.5)
    #timelines = {}
    #for i,strategy in enumerate(ylabels):
    #    timelines[strategy] = results[strategy]["cadence_hr"]
    fig = plt.figure(figsize=(20,8))
    ax = fig.add_subplot(111)
    marker_size = 50
    
    # Markers piechart
    # first define the ratios. Here one value because all the same
    ratio0 = 1./8   #  same for 8 filters
    ratio = ratio0
    # calculate the points of the first pie marker
    #
    # these are just the origin (0,0) +
    # some points on a circle cos,sin
    xy_list = []
    r1 = ratio0
    r2 = r1 + ratio0
    r3 = r2 + ratio0
    r4 = r3 + ratio0
    r5 = r4 + ratio0
    r6 = r5 + ratio0
    r7 = r6 + ratio0
    
    x = [0] + np.cos(np.linspace(0, 2*math.pi*r1, 10)).tolist()
    y = [0] + np.sin(np.linspace(0, 2*math.pi*r1, 10)).tolist()
    xy_list.append(list(zip(x, y)))
    
    x = [0] + np.cos(np.linspace(2*math.pi*r1, 2*math.pi*r2, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r1, 2*math.pi*r2, 10)).tolist()
    xy_list.append(list(zip(x, y)))

    x = [0] + np.cos(np.linspace(2*math.pi*r2, 2*math.pi*r3, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r2, 2*math.pi*r3, 10)).tolist()
    xy_list.append(list(zip(x, y)))

    x = [0] + np.cos(np.linspace(2*math.pi*r3, 2*math.pi*r4, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r3, 2*math.pi*r4, 10)).tolist()
    xy_list.append(list(zip(x, y)))

    x = [0] + np.cos(np.linspace(2*math.pi*r4, 2*math.pi*r5, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r4, 2*math.pi*r5, 10)).tolist()
    xy_list.append(list(zip(x, y)))

    x = [0] + np.cos(np.linspace(2*math.pi*r5, 2*math.pi*r6, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r5, 2*math.pi*r6, 10)).tolist()
    xy_list.append(list(zip(x, y)))    

    x = [0] + np.cos(np.linspace(2*math.pi*r6, 2*math.pi*r7, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r6, 2*math.pi*r7, 10)).tolist()
    xy_list.append(list(zip(x, y)))    

    x = [0] + np.cos(np.linspace(2*math.pi*r7, 2*math.pi, 10)).tolist()
    y = [0] + np.sin(np.linspace(2*math.pi*r7, 2*math.pi, 10)).tolist()
    xy_list.append(list(zip(x, y)))    
    
    for i in range(len(ylabels)):
        alpha = 1
        for xy, f in zip(xy_list, filters):
            timeline = [d for d in days if f in results[ylabels[i]][d][exp]]
            ax.plot(timeline, [(i*0.5)+0.5] * len(timeline), marker=(xy),
                    ms=marker_size, markerfacecolor=filt_color[f], markeredgecolor='k', linestyle='none',
                    label=f.replace("f", "F"))

    locsy, labelsy = plt.yticks(pos,ylabels)
    plt.setp(labelsy, fontsize = 16)
    
#    ax.axis('tight')
    ##ax.set_xlim([50,4e4])
    ##ax.set_ylim(ymin = -0.1, ymax = ilen*0.5+0.5)
    ax.grid(which='both', color='grey', linestyle=':')
    ax.set_xlabel(f"Time from merger (days)", fontsize=30)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(0.47, 1.2),
            ncol=8, fancybox=True, shadow=False, fontsize=18, framealpha=0.8, borderpad=1.5)
    
    #ax.set_xscale("log")
    ax.tick_params(labelsize=30, width=1, length=5)
    ax.set_ylim(ymin = 0.3, ymax = ilen*0.5+0.2)
    ax.invert_yaxis()
    # Fix the x the limits
    xlim = ax.get_xlim() # get existing x limits
    ax.set_xlim(xlim)
    # Force integers
    ax.xaxis.set_major_locator(MaxNLocator(integer=True)) 

    # PLOT horizontal line
    #ax.plot([xlim[0], xlim[1]], 2*[np.mean(np.arange(len(results))) - 0.25], color='grey')
    
    plt.subplots_adjust(left=0.2)

    #ax.legend(ncol=len(category_names), bbox_to_anchor=(0, 1),
    #          loc='lower left', fontsize='small')
    plt.savefig(f"gantt_chart_filters_{exp}_z{z}.pdf", bbox_inches='tight')

    return fig, ax


#########

lines_dict = {"Rubin": "--", "Roman":"-"}
maglim_dict_roman_1hr = {
"f062": 28.5, "f087": 28.2, "f106": 28.1, "f129": 28.0, "f158": 28.0, "f184": 27.4, "f213": 26.2, "f146": 28.4
}
maglim_dict_roman_55s = {
"f062": 25.5, "f087": 25.1, "f106": 25.1, "f129": 25.0, "f158": 24.9, "f184": 24.4, "f213": 23.7, "f146": 25.9
}

titles = [r"GW170817 best fit at 1Gpc, $\theta$=26 deg",
          r"GW170817 best fit at 1Gpc, $\theta$=90 deg",
          r"Kilonova low $Y_e$, $M_{ej}$ at 1Gpc, $\theta$=45 deg",
          r"Kilonova low $Y_e$ at 1Gpc, $\theta$=26 deg",
          r"Kilonova low $M_{ej}$ at 1Gpc, $\theta$=26 deg",
          r"Kilonova high $M_{ej}$ at 1Gpc, $\theta$=26 deg"]

out_filenames = ["figures_lc/plot_lc_gw170817_polar.pdf",
                 "figures_lc/plot_lc_gw170817_equatorial.pdf",
                 "figures_lc/plot_lc_low-all.pdf",
                 "figures_lc/plot_lc_low-ye.pdf",
                 "figures_lc/plot_lc_low-wind-ejecta.pdf",
                 "figures_lc/plot_lc_high-wind-ejecta.pdf"]

model_dict = {"figures_lc/plot_lc_gw170817_polar.pdf": "GW170817$_{pol}$",
              "figures_lc/plot_lc_gw170817_equatorial.pdf": "GW170817$_{eq}$",
              "figures_lc/plot_lc_low-all.pdf": "low $Y_{e}$, $M_{ej}$",
              "figures_lc/plot_lc_low-ye.pdf": "low $Y_{e}$",
              "figures_lc/plot_lc_high-wind-ejecta.pdf": "high $M_{ej}$",
              "figures_lc/plot_lc_low-wind-ejecta.pdf": "low $M_{ej}$"
}
"""
Bulla models: dyn has 3 values: mass, mean velocity, mean Ye
              dyn has 2 values: mass, mean velocity (fixed Ye)
"""
                
models = ["/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.005-0.15-0.20_wind0.010-0.05_theta25.84_z0.2.dat",
          "/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.005-0.15-0.20_wind0.010-0.05_theta90.00_z0.2.dat",
          "/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.001-0.15-0.15_wind0.010-0.05_theta45.57_z0.2.dat",
          "/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.005-0.15-0.15_wind0.010-0.05_theta25.84_z0.2.dat",
          "/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.001-0.15-0.20_wind0.010-0.05_theta25.84_z0.2.dat",
          "/Users/igor/data/POSSIS/POSSIS_lc_191022/nph1.0e+06_dyn0.010-0.15-0.20_wind0.090-0.05_theta25.84_z0.2.dat"]

def getDetections(t, days, filters):
    """
    Given a dictionary of limiting magnitudes, determine if the source is
    detected or not in the given filters

    Parameters
    ----------
    t astropy table
        light curve table
    days list
        list of days from the merger
    filters list
        list of filters to consider

    Return
    ------
    det_tot_shallow int
        total number of detections in 55s exp
    epochs_with_Ndet_shallow int
        number of epochs in which there are at least N detections in 55s exp
    det_tot_deep int
        total number of detections in 1hr exp
    result dict
        dictionary of filters for which there is a detection
    """
    epochs_with_Ndet_shallow = []
    det_tot_shallow = 0
    epochs_with_Ndet_deep = []
    det_tot_deep = 0
    # Results should be expressed in number of detections and filters per epoch
    result = {}
    # FIXME this is dumb
    for d in days:
        result[d] = {'55s': [], '1hr': []}
    # Iterate over the filters
    for filt in filters:
        # FIXME only Roman for now
        if filt[0] != "f":
            continue
        # Interpolate to reduce noise
        idx = [i for i in np.arange(len(t)) if np.isnan(t[filt][i]) == False]
        idx = [i for i in idx if not t[filt][i] == np.inf]
        #xnew = np.linspace(t["t[days]"][idx].min(), t["t[days]"][idx].max(), len(idx))
        f = interpolate.interp1d(t["t[days]"][idx],  t[filt][idx])    
        try:
            mags = f(days)
        except ValueError:
            mags = [float(f(d)) if (d <= np.max(t["t[days]"][idx])) else 99. for d in days]
        for i in np.arange(len(days)):
            if mags[i] <= maglim_dict_roman_1hr[filt]:
                det_tot_deep += 1
                result[days[i]]['1hr'].append(filt)
            if mags[i] <= maglim_dict_roman_55s[filt]:
                det_tot_shallow += 1
                result[days[i]]['55s'].append(filt)

    return det_tot_shallow, det_tot_deep, result


def getRawPotential(rate, area, time_window, maglim, M):
    """
    How many kilonovae shall we expect given an intrinsic
    BNS merger rate, an area, a time window, given our magnitude
    limit and the peak magnitude of the kilonova?

    Parameters
    ----------
    rate float
        in events/Gpc3/y
    area float
        in deg2
    time_window float
        in days
    maglim float
        magnitude limit (AB)
    Magpeak float
        peak absolute magnitude

    Returns
    -------
    n float
        number of expected events
    """
    # Find the redshift
    z = ac.z_at_value(cosmo.distmod, (maglim - M)*u.mag)
    vol = cosmo.comoving_volume(z).to("Gpc3")
    n = rate * (time_window / 365.) * vol.value * area/41253

    return n

def paper_getRawPotential():
    # Create table for the paper
    #BNS merger rate: 210^+240_-120 Gpc-3 y-1

    area = 19.
    time_window = 365 # 1 year
    rate = 210.
    rate_high = rate + 240
    rate_low = rate - 120
    Ms = [-14, -15, -16, -17]
    l = "Lim"
    for M in Ms:
        l += f" & M=${M}$"
    l += " \\" + "\\"
    print(l)
    print("\\hline")
    for maglim in [25., 28.]:
        l = f"{maglim} "
        for M in Ms:
            n = getRawPotential(rate, area, time_window, maglim, M)
            n_low = n - getRawPotential(rate_low, area, time_window, maglim, M)
            n_high = getRawPotential(rate_high, area, time_window, maglim, M) - n
            l += f" & ${'{:.2f}'.format(n)}^{{+{'{:.2f}'.format(n_high)}}}_{{-{'{:.2f}'.format(n_low)}}}$"
        l += " \\" + "\\"
        print(l)


def doPlot(t, out_filename, title, doShow=True, doSave=True, z=0.2):
        # Read the data file
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

        # Plot sensitivity
        if z < 0.1:
            ax.plot([0, 55], [20.5,20.5], color='grey', linewidth=3, linestyle=':', label="ZTF r limit 30s")
        else:
            ax.plot([0, 55], [22.2,22.2], color='grey', linewidth=3, linestyle=':', label="ZTF r limit 300s")
        ax.plot([0, 55], [24.03,24.03], color='k', linewidth=3, linestyle='--', label="LSST r limit 30s")
        ax.plot([0, 55], [25.0,25.0], color='orchid', linewidth=3, linestyle='-', label="F129 limit 55s")
        ax.plot([0, 55], [28.0,28.0], color='fuchsia', linewidth=3, linestyle='-', label="F129 limit 1hr")
        # And in the insets
        axins.plot([0, 55], [24.03,24.03], color='k', linewidth=3, linestyle='--', label="LSST r limit 30s")
        axins.plot([0, 55], [25.0,25.0], color='orchid', linewidth=3, linestyle='-', label="F129 limit 55s")

        # Set plot parameters
        plt.rcParams['xtick.labelsize']=11
        ax.legend(fontsize=11, loc='center left', bbox_to_anchor=(1, 0.5))
        if z == 0.05:
            ax.set_xlim([0.5, 27])
            ax.set_ylim([29, 19.5])
            # for the inset
            x1, x2, y1, y2 = 0.5, 7.5, 26, 19.5
        else:
            ax.set_xlim([0.5, 21])
            ax.set_ylim([29, 22.0])
            # for the inset
            x1, x2, y1, y2 = 0.5, 3.5, 26, 22.5
        ax.set_xlabel("Time from merger (days)", fontsize=13)
        ax.set_ylabel("AB Magnitude", fontsize=13)
        ax.set_title(title, fontsize=13)
        ax.tick_params(axis='both', labelsize=13)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

        # subregion of the original image
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        #axins.set_xticklabels([])
        axins.set_yticklabels([])
        #ax.indicate_inset_zoom(axins, edgecolor="black")
        if doSave is True:
            plt.savefig(out_filename, bbox_inches='tight')
        if doShow is True:
            plt.show()


def MultiDetectMetric(t, N=2, M=1, gap=4, filters=None):
    """
    get the fraction of sources detected at least N times in at least M filters
    given a certain cadence, here indicated as 'gap' in days
    """
    # random shift, determined by the gap itself
    days = np.linspace(1, 18-gap, int(np.round(18-gap)))
    days += np.random.random()*gap
    if filters is None:
        filters = t.colnames
    det_tot_shallow, det_tot_deep, result = getDetections(t, days, filters)
    # Filters with at least N detections
    filt_all = []
    for d in days:
        filt_all += result[d][exp]
    filters_det = [ff for ff in tab.colnames if len([x for x in filt_all if x==ff]) >= N]
    if len(filters_det) >= M:
        return 1
    else:
        return 0
    

def getFractionDetect(t, N=2, M=1, gap=4, n_iterations=100, filters=None):
    it = 0
    count = 0
    if filters is None:
        filters = t.colnames
    while it < n_iterations:
        it += 1
        count += MultiDetectMetric(t, N=N, M=M, gap=gap, filters=filters)

    return count/n_iterations


if __name__ == "__main__":
    getMetric = True
    Area = 19.
    N = 2 # number of detections
    M = 3 # in at least M filters
    gap = 6 # gap between observations
    Ntot = 100000 # Total number of injections
    distances_Mpc = np.array(list(np.linspace(100, 2000, 20)) + [3000, 4000, 5000, 6000])
    rate = 210.
    rate_high = rate + 240
    rate_low = rate - 120
    #selected_filters = ['f146']
    selected_filters = ['f062', 'f087', 'f106', 'f129', 'f158', 'f184', 'f213']
    # Raw potential table
    paper_getRawPotential()

    # redshift (available: 0.05, 0.02)
    #z = 0.05
    z = 0.2
    if z != 0.2:
        models = [x.replace("z0.2.dat", f"z{z}.dat") for x in models]
        titles = [x.replace("1Gpc", "230Mpc") for x in titles]
    out_filenames = [x.replace(".pdf", f"_z{z}.pdf") for x in out_filenames]

    if getMetric is False:
        print("\hline\hline")
        print("model & tot. & epochs & epochs &  epochs \
& filters  & filters & filters \\"+"\\")
        print("      & det. & $\geq$ 1 det & $\geq$ 2 det &  $\geq$ 3 det \
& 2 det  & 3 det & 4 det \\"+"\\")
        print("\hline")
    for exp in ['55s', '1hr']:
        if getMetric is True:
            print("EXPTIME:", exp)
            # Define distrubution between 100 and 1000 Mpc z = ac.z_at_value(cosmo.distmod, (maglim - M)*u.mag)
            z_list = [ac.z_at_value(cosmo.luminosity_distance, d) for d in distances_Mpc*u.Mpc]
            vols = [cosmo.comoving_volume(zz).to("Gpc3").value for zz in z_list]
            vols_norm = vols/np.sum(vols)
            # Assume what matters are all the params of filename but the distance
            for filename, out_filename, vol_norm in zip (models, out_filenames, vols_norm):
                frac_list = []
                frac0_list = []
                modelname = model_dict[out_filename.replace(f"_z{z}.pdf", ".pdf")]
                for d in distances_Mpc:
                    # Open the file for a certain model at a given distance
                    filename_d = filename.replace(f"z{z}", f"dMpc{int(d)}")
                    tab = ascii.read(filename_d)
                    # Run the metric
                    frac = getFractionDetect(tab, N=N, M=M, gap=gap,
                                             n_iterations=np.round(vol_norm*Ntot)+1,
                                             filters = selected_filters)
                    # Fraction for nightly cadence
                    frac0 = getFractionDetect(tab, N=N, M=M, gap=1.,
                                              n_iterations=np.round(vol_norm*Ntot)+1,
                                              filters = selected_filters)
                    frac_list.append(frac)
                    frac0_list.append(frac0)
                    #print(filename_d, frac)
                # Total numbers
                tot_kn = rate * np.sum([frac_list[i] * (vols[i] - vols[i-1]) for i in np.arange(len(vols)-1)+1]) * Area/41253
                tot_kn_high = rate_high * np.sum([frac_list[i] * (vols[i] - vols[i-1]) for i in np.arange(len(vols)-1)+1]) * Area/41253
                tot_kn_low = rate_low * np.sum([frac_list[i] * (vols[i] - vols[i-1]) for i in np.arange(len(vols)-1)+1]) * Area/41253
                tot_kn0 = rate * np.sum([frac0_list[i] * (vols[i] - vols[i-1]) for i in np.arange(len(vols)-1)+1]) * Area/41253
                print(f"{modelname} & ${'{:.1f}'.format(tot_kn)}^{{+{'{:.1f}'.format(tot_kn_high-tot_kn)}}}_{{-{'{:.1f}'.format(tot_kn-tot_kn_low)}}}$ & {'{:.1f}'.format(tot_kn/tot_kn0)}")
            raw_pot = rate * np.sum([1 * (vols[i] - vols[i-1]) for i in np.arange(len(vols)-1)+1]) * Area/41253
            print(f"N of kilonovae within {distances_Mpc[-1]} Mpc in an area of {Area}deg2: {rate*vols[-1] * Area/41253}/y")
            print("----")
        if getMetric is True:
            continue
        print(f"\multicolumn{{8}}{{c}}{{Exposure: {exp} }} \\"+"\\")
        print("\hline")
        results = {}
        for filename, out_filename, title in zip (models, out_filenames, titles):
            tab = ascii.read(filename)
            rate = 230
            area = 19
            time_window = 365
            maglim = 28.
            Magpeak = -16.
            #n = getRawPotential(rate, area, time_window, maglim, Magpeak)
            days = [3, 6, 9, 12, 15, 18]
            days = [2,4,6,8,10,12,14,16,18]
            det_tot_shallow, det_tot_deep, result = getDetections(tab, days, tab.colnames)
            if exp == '55s':
                det_tot = det_tot_shallow
            else:
                det_tot = det_tot_deep
            filt_all = []
            for d in days:
                filt_all += result[d][exp]
            modelname = model_dict[out_filename.replace(f"_z{z}.pdf", ".pdf")]
            print(modelname, "&",  det_tot, "&",
                  ",".join([str(k) for k in result.keys() if len(result[k][exp]) >= 1]), "&",
                  ",".join([str(k) for k in result.keys() if len(result[k][exp]) >= 2]), "&",
                  ",".join([str(k) for k in result.keys() if len(result[k][exp]) >= 3]), "&",
                  ",".join([ff for ff in tab.colnames if len([x for x in filt_all if x==ff]) == 2]), "&",
                  ",".join([ff for ff in tab.colnames if len([x for x in filt_all if x==ff]) == 3]), "&",
                  ",".join([ff for ff in tab.colnames if len([x for x in filt_all if x==ff]) == 4]), "\\"+"\\",
                  )
            results[modelname] = result


        if getMetric is False:
            figg, axg = chart([ff for ff in tab.colnames if ff[0] == 'f'], days, results, exp, z)
            plt.show()
            print("\hline")

    #for filename, out_filename, title in zip(models, out_filenames, titles):
    #    tab = ascii.read(filename)
    #    doPlot(tab, out_filename, title, doShow=True, doSave=True, z=z)

        
