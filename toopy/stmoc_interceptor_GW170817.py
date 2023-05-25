######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
import argparse
import time
import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time



# IMPORTS
import glob
import math
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.coordinates import Angle, SkyCoord
import astropy_healpix as ah
from mocpy import MOC
from mocpy import STMOC
from mocpy import WCS
from mocpy import WCS, STMOC, TimeMOC
from astroquery.vizier import Vizier
from astroquery.vizier import VizierClass
import pandas as pd
import healpy as hp
from astropy.time import TimeDelta
import numpy as np
import pandas as pd
import healpy as hp
import numpy as np
from urllib.request import urlopen
from astropy.utils.data import download_file


plt.switch_backend('agg')
######################################################################################
######################################################################################
# Args
######################################################################################
######################################################################################
start_fermipy = time.time()
parser = argparse.ArgumentParser(description="etc")
parser.add_argument("-string", "--string",
                    type=str,
                    help="string",
                    required=True)
args = parser.parse_args()

print(args.string)
outdir = './STMOC'
if not os.path.exists(outdir):
    os.mkdir(outdir)


base_dir=os.getcwd()
print(os.getcwd())
os.chdir(outdir)
print('######################################################################################')
print('######################################################################################')
print('################################################ GW ###############################')
print('######################################################################################')
df_gw = pd.read_csv('./AA_df_gw_STMOC_trial.csv', sep=',',  lineterminator='\n', names=None) 
print(df_gw)
stmoc_gw=STMOC.from_fits('./GW170817_&_Revision_1_STMOC.fits')
print("Time of the first gw observation: ", stmoc_gw.min_time.iso)
print("Time of the last observation: ", stmoc_gw.max_time.iso)


tmoc_gw = TimeMOC.from_time_ranges(
min_times = stmoc_gw.min_time, 
max_times = stmoc_gw.max_time,
)


moc_of_gw = stmoc_gw.query_by_time(tmoc_gw)


text_file = open('./random.txt', 'w')
text_file.write('BLa_gbm')
text_file.close()


print('######################################################################################')
print('######################################################################################')
print('################################################ GBM ###############################')
print('######################################################################################')
df_gbm = pd.read_csv('./AA_df_GRB170817A_STMOC_trial.csv', sep=',',  lineterminator='\n', names=None) 
print(df_gbm)
stmoc_gbm=STMOC.from_fits('./GRB170817A_FlightPos_STMOC.fits')
print("Time of the first gbm observation: ", stmoc_gbm.min_time.iso)
print("Time of the last observation: ", stmoc_gbm.max_time.iso)


tmoc_gbm = TimeMOC.from_time_ranges(
min_times = stmoc_gbm.min_time, 
max_times = stmoc_gbm.max_time,
)


moc_of_gbm = stmoc_gbm.query_by_time(tmoc_gbm)


#moc_gbm_MOC = stmoc_gbm.query_by_time(Time([[stmoc_gbm.min_time.iso, stmoc_gbm.max_time.iso]], format="iso", scale="tdb"))
text_file = open('./random.txt', 'w')
text_file.write('BLa_gbm')
text_file.close()
######################################################################################
######################################################################################
# Fermipy Analysis
######################################################################################
######################################################################################
print('######################################################################################')
print('######################################################################################')
print('################################################ intersect ###############################')
print('######################################################################################')
intersect_snt = stmoc_gw.intersection(stmoc_gbm)
print('Is the intersection empty? '+str(intersect_snt.is_empty()))

print("Time of the first observation: ", intersect_snt.min_time.iso)
print("Time of the last observation: ", intersect_snt.max_time.iso)





tmoc = TimeMOC.from_time_ranges(
min_times = intersect_snt.min_time, 
max_times = intersect_snt.max_time,
)


moc_of_intersect = intersect_snt.query_by_time(tmoc)

vizier = VizierClass(
row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
cat_4FGL, = vizier.get_catalogs('J/ApJS/247/33/4fgl')
moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
fig = plt.figure(111, figsize=(15, 10))
# Define a astropy WCS easily
with WCS(fig, 
        fov=330 * u.deg,
        center=SkyCoord(0,0, unit='deg', frame='galactic'),
        coordsys="galactic",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    moc_of_intersect.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="black", label='intersect')
    moc_of_intersect.border(ax=ax, wcs=wcs, alpha=0.5, color="black")

    moc_of_gw.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="red", label='GW170817')
    moc_of_gw.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    moc_of_gbm.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="grey", label='Fermi-GBM')
    moc_of_gbm.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    #moc_gbm_MOC.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green", label='GBM')
    #moc_gbm_MOC.border(ax=ax, wcs=wcs, alpha=0.5, color="black")

    moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
    moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
plt.xlabel('ra')
plt.ylabel('dec')
plt.legend(prop={'size': 20}, loc='upper right')
#plt.title('Query by time result')
plt.grid(color="black", linestyle="dotted")
outname = 'FOV_Galactic_GW&GRB_17.pdf'  
plt.savefig(outname)


###############
#Figure 2
###############
fig = plt.figure(figsize=(10, 10))
with WCS(fig, 
    fov=90 * u.deg,
    center=SkyCoord(df_gbm['RA_OBJ'][0], df_gbm['DEC_OBJ'][0], unit='deg', frame='icrs'),
    coordsys='icrs',
    rotation=Angle(0, u.degree),
    projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
moc_of_intersect.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="black", label='intersect')
moc_of_intersect.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
moc_of_gw.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="red", label='GW170817')
moc_of_gw.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
moc_of_gbm.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="grey", label='Fermi-GBM')
moc_of_gbm.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
ax.legend(prop={'size': 30}, loc='upper right')
plt.grid(color="black", linestyle="dotted")
plt.xlabel('RA', size=30)
plt.ylabel('DEC', size=30)
plt.xticks(color='black', fontsize=22)

ax.coords.grid(True, color='black')
#ax.coords[0].set_ticks(width=10)
ax.coords[0].set_ticklabel(size="xx-large")
#ax.coords[1].set_ticks(width=10)
ax.coords[1].set_ticklabel(size="xx-large")
outname = 'FOV_ROI_GW&GRB_17.pdf'
plt.savefig(outname)

print('######################################################################################')
print('######################################################################################')
print('#################### Archival ROI -- 4FGL ###############################')
print('######################################################################################')

print('######################################################################################')
print('######################################################################################')
print('Done and '+'Total Run-time for stmoc_interceptor.py is: '+str(time.time() - start_fermipy))
print('######################################################################################')
print('######################################################################################')


