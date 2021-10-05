#!/usr/bin/env python

import numpy as np
from matplotlib import gridspec
from matplotlib import font_manager
from matplotlib import pyplot as plt
from scipy import integrate
import sys
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)

if (sys.version_info[0] > 2):
    print("your python version is 3")
    from itertools import zip_longest as izip_longest
else:
    # Python 2 code in this block
    print("your python version is 2")
    from itertools import izip_longest


parser.add_argument('-emin', action='store', dest='emin',
                    help='Emin: Where the COHP becomes flat (default is -12)', type=float, required=False, default=-10)
parser.add_argument('-emax', action='store', dest='emax',
                    help='Emax: Where the COHP becomes flat (default is 8)', type=float, required=False, default=8)
parser.add_argument('-ymin', action='store', dest='ymin',
                    help='Do not set if not necessary', type=float, required=False)
parser.add_argument('-ymax', action='store', dest='ymax',
                    help='Do not set if not necessary', type=float, required=False)
parser.add_argument('--version', action='version', version='%(prog)s 1.00')
parser.add_argument('-print', action='store', dest='printed',
                    help='print datafile for each paoir (True/False---Default is True)', required=False, type=str, default="True")
parser.add_argument('-ratioscale', action='store', dest='ratioscale',
                    help='Scale the width by a defined value (default = 1)', required=False, type=float, default=1)

results = parser.parse_args()
Emin = results.emin
Emax = results.emax
prt = results.printed
ratio = results.ratioscale
#####################Plot Information############################
font_manager.findfont('Helvetica Light')
plt.rc('font', family='Helvetica Light')
plt.rc('font', serif='Helvetica Light', size=8)
# plt.rc('font', family='sans-serif')
plt.rcParams['axes.linewidth'] = 0.75
plt.rcParams['xtick.major.size'] = 4
plt.rcParams['xtick.minor.size'] = 2
plt.rcParams['xtick.major.width'] = 0.75
plt.rcParams['xtick.minor.width'] = 0.75
plt.rcParams['ytick.major.size'] = 4
plt.rcParams['ytick.minor.size'] = 2
plt.rcParams['ytick.major.width'] = 0.75
plt.rcParams['ytick.minor.width'] = 0.75
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.direction'] = 'in'
###################################################################

inputcohp = "COHPCAR.lobster"
collect = []
with open(inputcohp, 'r') as f:
    f.readline()  # just skip the first  line
    numcohp = (f.readline().split()[0])  # just skip the second line
    f.readline()  # just skip the third  line
    for j in range(int(numcohp)-1):
        collect.append(f.readline().replace('->', ' ').replace(')',
                                                               ' ').replace('(', ' ').split(':')[1:])
    data = list(izip_longest(*[x.split() for x in f]))
    data = np.array(data)
    datalife = data.astype(float)  # np.float is deprecated
    datalife = datalife

datafile = np.array(datalife)

pairs_tot = []
for j in range(len(collect)):
    pairs_tot.append(collect[j][0].split())

X = datafile[0]

# numcohp=2
fop = open('LOBSTER_INTEGRATION.DAT', 'w')
fop.write('#  Pairs    Dist [A]  BondOrder\n')
fop.write('# -----------------------------\n')

tot_col_in_cohp = (int(numcohp)*2)
for k in range(3, tot_col_in_cohp+1, 2):
    ########################################################
    fig = plt.figure(figsize=(3.23606797749978969640/ratio, 2))
    gs = gridspec.GridSpec(1, 1, width_ratios=[3])

    ax0 = fig.add_subplot(111)
    ax1 = fig.add_subplot(gs[0])
    ax0.spines['top'].set_color('none')
    ax0.spines['bottom'].set_color('none')
    ax0.spines['left'].set_color('none')
    ax0.spines['right'].set_color('none')
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax1.axhline(y=0, ls='-.', lw=0.5, color='k')
    ax1.axvline(x=0, ls='--', lw=0.5, color='gray')

    cohpcolor = '#e76f51'
    cohpfilll = '#e76f51'
    intcohp = '#8ecae6'
    intcohp_mine = '#023047'
    ########################################################
    X = datafile[0]
    Y = datafile[k]
    Z = datafile[k+1]
    SUMM = np.sum(np.array(Y)[(X > Emin) & (X <= 0)])*(X[1]-X[0])
    Z_int = integrate.cumtrapz(Y, X, initial=0)
    Z_intmin = np.amin(np.array(-Z_int)[(X > Emin) & (X < 0)])
    SUMM = -round(SUMM, 1)
    # Create Index to finding the element
    index = int((k-3)/2)
    dist = str(round(float(pairs_tot[index][2]), 2))
    name = str(pairs_tot[index][0])+'--'+str(pairs_tot[index][1])

    if prt == 'True':
        f = open(name+".dat", "w")
        f.write(
            '#The bond order (-intg) of between {} is about {}\n'.format(name, SUMM))
        f.write('#     Ene           COHP      Int(COHP)  Int(My integ.)\n')
        f.write('#------------------------------------------------------\n')
        for j in range(len(X)):
            f.write('%12.5f %12.5f %12.5f %12.5f\n' %
                    (X[j], Y[j], Z[j], Z_int[j]))
    fop.write('%10s %8.4f  %6.2f\n' % (str(name), float(dist), float(SUMM)))

    Zmin = np.amin(np.array(-Z)[(X > Emin) & (X < 0)])
    ax1.plot(X, -Y, color=cohpcolor, lw=0.25)
    ax1.fill_between(X, -Y, facecolor=cohpfilll, interpolate=True,
                     alpha=1, label=name+' ('+dist+r' $\AA$'+')')
#    ax1.plot(X, -Z-Zmin, color=intcohp, lw=0.5, label='ICOHP')
    ax1.plot(X, -Z_int-Z_intmin, color=intcohp_mine,
             lw=.5, ls='--', label='ICOHP='+str(SUMM))
#    ax1.plot([0, 0], [0, 0], lw=0, label='ICOHP='+str(SUMM))

    handles, labels = ax1.get_legend_handles_labels()
    handles = [handles[0], handles[1]]
    labels = [labels[0], labels[1]]
    ax1.legend(handles, labels, fontsize=6)

#    f.write("%s %18.2f\n" % (str(name), float(dist)))
#            fout.write("%8d %18.8f %18.8f\n" % (i+1, kdist2[i], m[i]))
#    ax1.legend(fontsize=6)
    # Emin=-5
    #####PLOT SETTING AXIS, Legend,etc.##############################
    YZ = np.concatenate((Y, Z), axis=None)
    XX = np.concatenate((X, X), axis=None)
    SelectedRange = np.array(YZ)[(XX > Emin) & (XX <= Emax)]

    if results.ymin:
        ymin = results.ymin
    else:
        ymin = -np.amax(SelectedRange)

    if results.ymax:
        ymax = results.ymax
    else:
        ymax = -np.amin(SelectedRange)

    ax1.set_xlim(Emin, Emax)
    ax1.set_ylim(ymin*1.05, ymax*1.05)
    ax1.set_ylabel(r'$-$COHP')
    ax1.set_xlabel(r'$E-E_{\rm F}$ (eV)')
    #################################################################
    plt.savefig(name+'.pdf', bbox_inches='tight', pad_inches=0.02)
#   plt.show()
    fig.clf()