######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
#basics
import matplotlib.pyplot as plt
import healpy as hp
import pandas as pd
import numpy as np
import math
import time
import os
from astropy.table import Table
from mpl_toolkits.mplot3d import Axes3D
#astropy
import astropy.units as u
import astropy_healpix as ah
from astropy.io import fits
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.utils.data import download_file
from astropy.coordinates import Angle, SkyCoord
#mocpy
from mocpy import MOC
from mocpy import STMOC
from mocpy import WCS
#vizier
from astroquery.vizier import VizierClass
#helper
from library.helper import observability_gbm
from library.helper import observability_gbm_gladePlus
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(event, trans_Num, vol_percent, rank, too_span, t_res, zenith, catalog):
        ###############
        #Event
        ###############
        filename = download_file(event, cache=True)
        hdul1 = fits.open(event)
        d_gbm = {'file': str(filename), 'URL': str(event),
               'trans_Num': str(trans_Num),
               'RA_OBJ': hdul1[0].header['RA_OBJ'],
               'DEC_OBJ':hdul1[0].header['DEC_OBJ'],
               'ERR_RAD':hdul1[0].header['ERR_RAD'],
               'DATE-OBS': hdul1[0].header['DATE-OBS'],
                'DATE-END': hdul1[0].header['DATE-END']}
        df_gbm = pd.DataFrame(data=d_gbm,index=[0])
        START=hdul1[0].header['DATE-OBS']
        print('#########################')
        print('#########################')
        print('This is START:'+str(START))
        print(trans_Num)
        print('#########################')
        print('#########################')
        skymap_event, header = hp.read_map(filename, h=True, verbose=False)
        quantile = vol_percent # to get x% contour
        header = dict(header)
        NSIDE = header['NSIDE']
        argsort = np.argsort(-skymap_event)
        cum_skymap_event = np.cumsum(sorted(skymap_event,reverse=True))
        cont_ind = argsort[cum_skymap_event < quantile]
        contour = np.array([1. if pix in cont_ind else 0. for pix in range(len(skymap_event))])
        max_pix = np.argmax(contour)
        wh_c=np.where(contour == 1)[0]
        wh2_c=np.argwhere(contour == 1)
        dec, ra = hp.pix2ang(nside=NSIDE, ipix=[max_pix])
        dec = math.pi/2. - dec
        theta, phi = hp.pix2ang(NSIDE, wh_c)
        ra_Cascade_evt = np.rad2deg(phi)
        dec_Cascade_evt = np.rad2deg(0.5 * np.pi - theta)
        skycoord_Cascade_evt = SkyCoord(ra_Cascade_evt, dec_Cascade_evt, unit="deg", frame="icrs")
        ###############
        #Rankings
        ###############
        if rank == 'Xmatch':
            outdir = './GBM_Alert/Xmatch/'+str(too_span)+'_ToO_TRes'+str(t_res)+str('hrs')+'_trigger_'+str(START)+'_&_'+str(trans_Num)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        ###############
        #4FGL
        ###############
        vizier = VizierClass(
        row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
        cat_4FGL, = vizier.get_catalogs(catalog)
        ###############
        #MOC
        ###############
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        moc_GBM_evt = MOC.from_skycoords(skycoord_Cascade_evt, max_norder=7)
        df_gbm['MOC']=moc_GBM_evt
        
        skycoord_evt = SkyCoord(hdul1[0].header['RA_OBJ'], hdul1[0].header['DEC_OBJ'], unit="deg", frame='icrs')
        SRC_ERROR50=hdul1[0].header['ERR_RAD']
        moc_evt_cone = MOC.from_cone(
        lon=skycoord_evt.ra,
        lat=skycoord_evt.dec,
        radius=Angle(SRC_ERROR50, u.deg),
        max_depth=10)
        df_gbm['MOC_cone']=moc_evt_cone
    
        ###############
        #Figure 1
        ###############
        fig = plt.figure(figsize=(10,10))
        with WCS(fig, 
                fov=330 * u.deg,
                center=SkyCoord(0,0, unit='deg', frame='galactic'),
                coordsys='galactic',
                rotation=Angle(0, u.degree),
                projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
            moc_GBM_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='HEALPIX (Contour: '+str(vol_percent*100)+' %)')
            moc_GBM_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
            moc_evt_cone.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='CONE (Error: '+str(SRC_ERROR50)+' deg)')
            moc_evt_cone.border(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="black")
            moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
            moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 20}, loc='upper right')
        plt.grid(color="black", linestyle="dotted")
        outname = 'FOV_Galactic.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        ###############
        #Figure 2
        ###############
        fig = plt.figure(figsize=(10, 10))
        with WCS(fig, 
            fov=30 * u.deg,
            center=SkyCoord(skycoord_evt.ra,skycoord_evt.dec, unit='deg', frame='icrs'),
            coordsys='icrs',
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc_GBM_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='HEALPIX (Contour: '+str(vol_percent*100)+' %)')
        moc_GBM_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
        moc_evt_cone.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='CONE (Error: '+str(SRC_ERROR50)+' deg)')
        moc_evt_cone.border(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="black")
        moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
        moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 20}, loc='upper right')
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
        ###################
        #XMatch with Glade2
        ###################
        cat_vizer = moc_GBM_evt.query_vizier_table('VII/281/glade2')
        data=pd.DataFrame({'RA': cat_vizer['_RAJ2000'], 'DEC':cat_vizer['_DEJ2000'],'dist':cat_vizer['Dist'],'Bmag':cat_vizer['Bmag'],'HyperLEDA':cat_vizer['HyperLEDA']})
        # ADD ME!!! #msk1=data[['dist']]<=distmax
        # ADD ME!!! #msk2=data[['dist']]>=distmin
        filter_good_ones =  (data.dist > 0) & \
                              (data.dist !='NaN') & \
                              (data.Bmag !='NaN') & \
                              (data.Bmag !='null') & \
                              (data.Bmag > 0)
        crossmatched_glade2_cat = data[filter_good_ones]
        outname = 'Xmatched_list_with_Glade2.csv'
        fullname = os.path.join(outdir, outname)    
        crossmatched_glade2_cat.to_csv(fullname, sep="\t", index = False, header=True)
        return crossmatched_glade2_cat, hdul1, outdir

    def Xmatched_raw_to_3Dplot_glade2(crossmatched_cat, hdul1, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat=crossmatched_cat.sort_values(by='Bmag', ascending=False)
        crossmatched_cat_top3=crossmatched_cat.head(3)
        print(crossmatched_cat_top3)
        print(hdul1[0].header)
        #data=Table.read('test_data.fits')
        #min_red=min(data['redshift'])
        '''fig = plt.figure(figsize=(16,14))
        ax = Axes3D(fig)

        ax = fig.gca(projection='3d')
        ax.view_init(10,90)

        z = crossmatched_cat.dist.values
        x = crossmatched_cat.RA.values
        y = crossmatched_cat.DEC.values
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
        z = crossmatched_cat.dist.values
        x = crossmatched_cat.RA.values
        y = crossmatched_cat.DEC.values
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
        hist, xedges, yedges = np.histogram2d(crossmatched_cat.RA.values, crossmatched_cat.DEC.values, bins=N_bins)
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

   
    def Xmatched_top10_BMag_to_obslist_glade2(event, observatory, crossmatched_cat, zenith, moon_sep, hdul1, too_span, time_resolution, outdir):
        ###############
        #Event
        ###############
        filename = download_file(event, cache=True)
        hdul1_n = fits.open(event)
        print(hdul1_n)
        print('This is zenith:'+str(zenith))
        print('This is time_resolution:'+str(time_resolution))
        crossmatched_cat=crossmatched_cat.sort_values(by='Bmag', ascending=False)
        crossmatched_cat=crossmatched_cat.head(10)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_gbm.merged_def2.doit(observatory, crossmatched_cat, zenith, moon_sep, hdul1, too_span, time_resolution, outdir)
        ###############
        #Pandas
        ###############
        listed_obs=[]
        for i in range(0, len(c_fin)):
            dict = {'RA': crossmatched_cat.RA.values[i],
            'DEC': crossmatched_cat.DEC.values[i],
            'dist': crossmatched_cat.dist.values[i],
            'Bmag': crossmatched_cat.Bmag.values[i],
            'HyperLEDA': crossmatched_cat.HyperLEDA.values[i],
            'Observable?': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Observatory': observatory,
            'Timeslot':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Airmass': airmass,
            'Altitude': altitude,
            'Zenith': zenith}
            #'offset': crossmatched_cat.offset.values[i]}     
            observability_df=pd.DataFrame(dict, index = [item for item in c_fin[i]])
            listed_obs.append(observability_df)
        listed_obs=pd.concat(listed_obs)
        fin_df=listed_obs.loc[True]
        outname = 'Glade2_Observability_@'+str(observatory)+'.csv'
        fullname = os.path.join(outdir, outname)    
        fin_df.to_csv(fullname, sep="\t", index = False, header=True)

        