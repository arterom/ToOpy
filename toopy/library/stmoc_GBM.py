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
        #if rank == 'Xmatch':
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
        #STMOC
        ###############
        outdir_stmoc = './STMOC'
        if not os.path.exists(outdir_stmoc):
            os.mkdir(outdir_stmoc)
        string_list_gbm=[]
        for i in range(0,len(df_gbm)):
            string=str(df_gbm['DATE-OBS'].values[i])
            string_list_gbm.append(string)
        times_gbm = Time(string_list_gbm, format='isot', scale='utc')
        dt_iso_gbm = TimeDelta(200, format='sec')
        t_gbm_start=times_gbm-dt_iso_gbm
        t_gbm_end=times_gbm+dt_iso_gbm
        df_gbm['STMOC_gbm_start']=t_gbm_start
        df_gbm['STMOC_gbm_stop']=t_gbm_end

        output_path='./STMOC/AA_df_gbm_STMOC_trial.csv'
        df_gbm.to_csv(output_path, mode='a', header=not os.path.exists(output_path))

        stmoc_gbm = STMOC.from_spatial_coverages(t_gbm_start, t_gbm_end, df_gbm['MOC_cone'])
        print("Time of the first observation: ", stmoc_gbm.min_time.iso)
        print("Time of the last observation: ", stmoc_gbm.max_time.iso)

        output_path='./STMOC/GBM_TransNum_'+str(trans_Num)+'_STMOC.fits'
        stmoc_gbm.write(output_path,format='fits', overwrite=True)

        #Stacked
        #stacked_STMOC=STMOC.from_fits('./GBM_Alert/stacked_STMOC.fits')
        #restacked_STMOC=stmoc_gbm.union(stacked_STMOC)
        #output_path='./GBM_Alert/stacked_STMOC.fits'
        #restacked_STMOC.write(output_path,format='fits', overwrite=True)

        if os.path.isfile('./STMOC/AA_stacked_STMOC.fits'):
            print ("File does exist")
            stacked_STMOC=STMOC.from_fits('./STMOC/AA_stacked_STMOC.fits')
            restacked_STMOC=stmoc_gbm.union(stacked_STMOC)
            output_path='./STMOC/AA_stacked_STMOC.fits'
            restacked_STMOC.write(output_path,format='fits', overwrite=True)
        else:
            print ("File does not exist")
            output_path='./STMOC/AA_stacked_STMOC.fits'
            stmoc_gbm.write(output_path,format='fits')
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

   
   
