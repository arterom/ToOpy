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
print('################################################ SWIFT ###############################')
print('######################################################################################')
df_swift = pd.read_csv('./AA_df_swift_STMOC_trial.csv', sep=',',  lineterminator='\n', names=None) 
print(df_swift)
stmoc_swift=STMOC.from_fits('./Swift_TrigID_1142847_STMOC.fits')
print("Time of the first swift observation: ", stmoc_swift.min_time.iso)
print("Time of the last observation: ", stmoc_swift.max_time.iso)

tmoc_swift = TimeMOC.from_time_ranges(
min_times = stmoc_swift.min_time, 
max_times = stmoc_swift.max_time,
)


moc_of_swift = stmoc_swift.query_by_time(tmoc_swift)

#moc_swift_MOC = stmoc_swift.query_by_time(Time([[stmoc_swift.min_time.iso, stmoc_swift.max_time.iso]], format="iso", scale="tdb"))
text_file = open('./random.txt', 'w')
text_file.write('BLa_Swift')
text_file.close()

print('######################################################################################')
print('######################################################################################')
print('################################################ GBM ###############################')
print('######################################################################################')
df_gbm = pd.read_csv('./AA_df_gbm_STMOC_trial.csv', sep=',',  lineterminator='\n', names=None) 
print(df_gbm)
stmoc_gbm=STMOC.from_fits('./GBM_TransNum_691590290_STMOC.fits')
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
intersect_snt = stmoc_swift.intersection(stmoc_gbm)
print(str(intersect_snt.is_empty()))

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
        center=SkyCoord(0, 0, unit='deg', frame='galactic'),
        coordsys="galactic",
        rotation=Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    moc_of_intersect.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="blue", label='intersect')
    moc_of_intersect.border(ax=ax, wcs=wcs, alpha=0.5, color="black")

    moc_of_swift.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="red", label='Swift')
    moc_of_swift.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    moc_of_gbm.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='GBM')
    moc_of_gbm.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    #moc_gbm_MOC.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green", label='GBM')
    #moc_gbm_MOC.border(ax=ax, wcs=wcs, alpha=0.5, color="black")

    moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
    moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
plt.xlabel('ra')
plt.ylabel('dec')
plt.legend()
#plt.title('Query by time result')
plt.grid(color="black", linestyle="dotted")
outname = 'FOV_Galactic.pdf'  
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


