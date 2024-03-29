######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
#basics
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import os
#astropy
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord
#mocpy
from mocpy import MOC
from mocpy import STMOC
from mocpy import WCS
#vizier
from astroquery.vizier import VizierClass
from astroquery.ned import Ned
#helper
from library.helper import observability_swift

from astropy.time import Time
from astropy.time import TimeDelta
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(too_span, t_res, catalog, ra, dec, err, night, TrigID):
        ###############
        #Event
        ###############
        observing_night=Time(night, format='mjd', scale='utc').isot


        #TESTING!!!!!!!!!!!!!!    
        #start_time_mjd = Time(night, format='mjd', scale='utc')
        #dt_sod = TimeDelta(46827.12, format='sec')
        #start_time_mjd_actual = start_time_mjd + dt_sod
        #observing_night_actual=Time(start_time_mjd_actual, format='mjd', scale='utc').isot
        #TESTING!!!!!!!!!!!!!!   
        d_swift = {
               'TrigID': str(TrigID),
               'RA_OBJ': ra,
               'DEC_OBJ':dec,
               'ERR_RAD':err,
               'DATE-OBS': observing_night,
               'MJD-OBS': night}
               #'start_time_mjd_actual': start_time_mjd_actual,
               #'observing_night_actual': observing_night_actual}
        df_swift = pd.DataFrame(data=d_swift,index=[0])
        print('#############################################')
        print('#############################################')
        print(df_swift)
        print('#############################################')
        print('###############')

        print('Manual input of alert date: '+str(observing_night))
        outdir = './Swift_BAT_Alert/Xmatch/'+str(too_span)+'_ToO_TRes'+str(t_res)+str('hrs')+'_&_'+str(observing_night)+'_&_TrigID_'+str(TrigID)
        if not os.path.exists(outdir):
            os.mkdir(outdir)









        ###############
        #4FGL-catalog
        ###############
        vizier = VizierClass(
        row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
        cat_4FGL, = vizier.get_catalogs(catalog)
        ###############
        #MOC
        ###############
        skycoord_evt = SkyCoord(ra, dec, unit='deg', frame='icrs')
        moc_Swift_BAT = MOC.from_cone(
        lon=skycoord_evt.ra,
        lat=skycoord_evt.dec,
        radius=Angle(err, u.deg),
        max_depth=10
        )
        df_swift['MOC']=moc_Swift_BAT
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        ###############
        #Glade2
        ###############
        cat_vizer = moc_Swift_BAT.query_vizier_table('VII/281/glade2')
        data=pd.DataFrame({'RA': cat_vizer['_RAJ2000'], 'DEC':cat_vizer['_DEJ2000'],'dist':cat_vizer['Dist'],'Bmag':cat_vizer['Bmag'],'HyperLEDA':cat_vizer['HyperLEDA']})
        # ADD ME!!! #msk1=data[['dist']]<=distmax
        # ADD ME!!! #msk2=data[['dist']]>=distmin
        filter_good_ones =  (data.dist > 0) & \
                              (data.dist !='NaN') & \
                              (data.Bmag !='NaN') & \
                              (data.Bmag !='null') & \
                              (data.Bmag > 0)
        crossmatched_cat_glade2 = data[filter_good_ones]
        outname = 'Xmatched_list_Glade2.csv'
        fullname = os.path.join(outdir, outname)    
        crossmatched_cat_glade2.to_csv(fullname, sep="\t", index = False, header=True)

        print('################################################################################################')
        print('##################### crossmatched_cat_glade2 #########################################################')
        print('################################################################################################')
        print(crossmatched_cat_glade2)
        print('################################################################################################')
        print('################################################################################################')
        
        ###############
        #Figure 1
        ###############
        fig = plt.figure(figsize=(10, 10))
        with WCS(fig, 
            fov=330 * u.deg,
            center=SkyCoord(0,0, unit='deg', frame='galactic'),
            coordsys='galactic',
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc_Swift_BAT.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='SwiftBAT Event')
        moc_Swift_BAT.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
        moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
        moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 30}, loc='upper right')
        plt.grid(color="black", linestyle="dotted")
        outname = 'FOV_Galactic.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        ###############
        #Figure 2
        ###############
        fig = plt.figure(figsize=(10, 10))
        with WCS(fig, 
            fov=15 * u.deg,
            center=SkyCoord(skycoord_evt.ra,skycoord_evt.dec, unit='deg', frame='icrs'),
            coordsys='icrs',
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc_Swift_BAT.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='SwiftBAT Event')
        moc_Swift_BAT.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
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
        outname = 'FOV_ROI.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        return crossmatched_cat_glade2, outdir

    def Xmatched_raw_to_3Dplot(crossmatched_cat_glade2, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat_glade2=crossmatched_cat_glade2.sort_values(by='Bmag', ascending=False)
        crossmatched_cat_glade2_top3=crossmatched_cat_glade2.head(3)
        print(crossmatched_cat_glade2_top3)
        #data=Table.read('test_data.fits')
        #min_red=min(data['redshift'])
        '''fig = plt.figure(figsize=(16,14))
        ax = Axes3D(fig)

        ax = fig.gca(projection='3d')
        ax.view_init(10,90)

        z = crossmatched_cat_glade2.dist.values
        x = crossmatched_cat_glade2.RA.values
        y = crossmatched_cat_glade2.DEC.values
        min_dist=min(z)

        #y=list(data['dec'])
        ax.scatter(x, y, z,'ko', c=y, cmap = 'Greys')


        m=ax.plot(x, y, 'ro', markersize=.5, color='r', zdir='z', zs=min_dist)

        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('Distance')
        outname = 'Glade2_Dist.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)'''
        ##########################################################################
        ##########################################################################
        ##########################################################################
        ###############
        #Figure 2
        ###############      
        fig = plt.figure(figsize=(20,10))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.view_init(20,-45)
        z = crossmatched_cat_glade2.dist.values
        x = crossmatched_cat_glade2.RA.values
        y = crossmatched_cat_glade2.DEC.values
        min_dist=min(z)
        ax.scatter(x, y, z, color='black', alpha=.9)
        m=ax.plot(x, y, 'ro', markersize=.5, color='r', zdir='z', zs=min_dist)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('Distance')
        ##################################################################################
        ##################################################################################
        N_bins = 30
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax.view_init(20,-45)
        hist, xedges, yedges = np.histogram2d(crossmatched_cat_glade2.RA.values, crossmatched_cat_glade2.DEC.values, bins=N_bins)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)
        xpos = xpos.flatten('F')
        ypos = ypos.flatten('F')
        zpos = np.zeros_like(xpos)
        dx = 0.5 * np.ones_like(zpos)
        dy = dx.copy()
        dz = hist.flatten()
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='red', zsort='average')
        outname = 'Glade2_Dist_full_catalog.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        return outdir


    def Xmatched_to_obslist(observatory, crossmatched_cat_glade2, zenith, moon_sep, too_span, time_resolution, night, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat_glade2=crossmatched_cat_glade2
        print('########################')
        print('########################')
        print(crossmatched_cat_glade2)
        print('########################')
        print('########################')
        print('This is zenith:'+str(zenith))
        print('This is time_resolution:'+str(time_resolution))
        crossmatched_cat_glade2=crossmatched_cat_glade2.sort_values(by='Bmag', ascending=False)
        crossmatched_cat_glade2=crossmatched_cat_glade2.head(10)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_swift.merged_def2.doit(observatory, crossmatched_cat_glade2, zenith, moon_sep, night, too_span, time_resolution, outdir)
        ###############
        #Pandas
        ###############
        listed_obs=[]
        for i in range(0, len(c_fin)):
            dict = {'RA': crossmatched_cat_glade2.RA.values[i],
            'DEC': crossmatched_cat_glade2.DEC.values[i],
            'dist': crossmatched_cat_glade2.dist.values[i],
            'Bmag': crossmatched_cat_glade2.Bmag.values[i],
            'HyperLEDA': crossmatched_cat_glade2.HyperLEDA.values[i],
            'Observable': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Observatory': observatory,
            'Timeslot':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Airmass': airmass,
            'Altitude': altitude,
            'Zenith': zenith}
            observability_df=pd.DataFrame(dict, index = [item for item in c_fin[i]])
            listed_obs.append(observability_df)
        listed_obs=pd.concat(listed_obs)
        fin_df=listed_obs.loc[True]
        outname = 'Observability_@'+str(observatory)+'.csv'
        fullname = os.path.join(outdir, outname)    
        fin_df.to_csv(fullname, sep="\t", index = False, header=True)
