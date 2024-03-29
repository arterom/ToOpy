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
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(event, graceid, rev, vol_percent, rank,too_span, t_res, zenith, catalog):
        ###############
        #Event
        ###############
        start_download = time.time()
        filename = download_file(event, cache=True)
        hdul1 = fits.open(event)
        hdul1.info()
        hdul1[1].columns
        data = hdul1[1].data
        header = hdul1[1].header
        ORDERING = header['ORDERING']
        header = dict(header)
        print(header)
        print('######################################################################################')
        print('######################################################################################')
        print('######################################################################################')
        print('######################################################################################')
        MJD_OBS = header['MJD-OBS']
        DATE_OBS = header['DATE-OBS']
        DISTMEAN = header['DISTMEAN']
        DISTSTD = header['DISTSTD']
        DISTMIN = DISTMEAN-DISTSTD
        DISTMAX = DISTMEAN+DISTSTD
        out_directory_date = header['DATE-OBS'][:11]

        d_gw = {'file': str(filename), 'URL': str(event),
               'graceid': str(graceid),
               'Revision': str(rev),
               'MJD_OBS': MJD_OBS, 
               'DATE_OBS': DATE_OBS}
        df_gw = pd.DataFrame(data=d_gw,index=[0])
        print(df_gw)
        print('download done and '+'Total Run-time is: '+str(time.time() - start_download))
        start_moc_cat = time.time()
        ###############
        #Rankings
        ###############
        outdir = './GW_Alert/Xmatch/'+str(too_span)+'_ToO_TRes'+str(t_res)+str('hrs')+'_&_'+str(out_directory_date)+'_&_'+str(graceid)+'_&_Revision_'+str(rev)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        ###############
        #4FGL
        ###############
        vizier = VizierClass(
        row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
        cat_4FGL, = vizier.get_catalogs(catalog)
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        print('moc_cat done and '+'Total Run-time is: '+str(time.time() - start_moc_cat))
        start_moc_evt = time.time()
        if ORDERING == 'NESTED':
            NSIDE = header['NSIDE']
            skymap_event, header = hp.read_map(filename, h=True, verbose=False)
            
            print('######################################################################################')
            print('######################################################################################')
            print('######################################################################################')
            print('######################################################################################')
            print('This step takes some time......')
            quantile = vol_percent # to get x% contour
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
            ra_evt = np.rad2deg(phi)
            dec_evt = np.rad2deg(0.5 * np.pi - theta)
            skycoord_evt = SkyCoord(ra_evt, dec_evt, unit="deg", frame="icrs")
            moc_90_nside = MOC.from_skycoords(skycoord_evt, max_norder=7)
            df_gw['MOC']=moc_90_nside

            print('moc_evt done and '+'Total Run-time is: '+str(time.time() - start_moc_evt))
            start_stmoc_evt = time.time()
            ###############
            #STMOC
            ###############
            outdir_stmoc = './STMOC'
            if not os.path.exists(outdir_stmoc):
                os.mkdir(outdir_stmoc)
            string_list_gw=[]
            for i in range(0,len(df_gw)):
                string=str(df_gw['DATE_OBS'].values[i])
                string_list_gw.append(string)
            times_gw = Time(string_list_gw, format='isot', scale='utc')
            dt_iso_gw = TimeDelta(200, format='sec')
            t_GW_start=times_gw-dt_iso_gw
            t_GW_end=times_gw+dt_iso_gw
            df_gw['STMOC_GW_start']=t_GW_start
            df_gw['STMOC_GW_stop']=t_GW_end

            output_path='./STMOC/AA_df_gw_STMOC_trial.csv'
            df_gw.to_csv(output_path, mode='a', header=not os.path.exists(output_path))

            stmoc_gw = STMOC.from_spatial_coverages(t_GW_start, t_GW_end, df_gw['MOC'])
            print("Time of the first observation: ", stmoc_gw.min_time.iso)
            print("Time of the last observation: ", stmoc_gw.max_time.iso)

            output_path='./STMOC/'+str(graceid)+'_&_Revision_'+str(rev)+'_STMOC.fits'
            stmoc_gw.write(output_path,format='fits', overwrite=True)


            if os.path.isfile('./STMOC/AA_stacked_STMOC.fits'):
                print ("File does exist")
                stacked_STMOC=STMOC.from_fits('./STMOC/AA_stacked_STMOC.fits')
                restacked_STMOC=stmoc_gw.union(stacked_STMOC)
                output_path='./STMOC/AA_stacked_STMOC.fits'
                restacked_STMOC.write(output_path,format='fits', overwrite=True)
            else:
                print ("File does not exist")
                output_path='./STMOC/AA_stacked_STMOC.fits'
                stmoc_gw.write(output_path,format='fits')
            print('stmoc_evt done and '+'Total Run-time is: '+str(time.time() - start_stmoc_evt))
            ###############
            #Figure NESTED
            ###############
            fig = plt.figure(figsize=(10,10))
            with WCS(fig, 
                    fov=330 * u.deg,
                    center=SkyCoord(0,0, unit='deg', frame='galactic'),
                    coordsys='galactic',
                    rotation=Angle(0, u.degree),
                    projection="AIT") as wcs:
                ax = fig.add_subplot(1, 1, 1, projection=wcs)
                moc_90_nside.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='HEALPIX (Contour: '+str(vol_percent*100)+' %)')
                moc_90_nside.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
                moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
                moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            ax.legend(prop={'size': 10})
            plt.grid(color="black", linestyle="dotted")
            outname = 'FOV_Galactic_nside.pdf'
            fullname = os.path.join(outdir, outname)    
            plt.savefig(fullname)
            cat_vizer = moc_90_nside.query_vizier_table('VII/281/glade2')
            data=pd.DataFrame({'RA': cat_vizer['_RAJ2000'], 'DEC':cat_vizer['_DEJ2000'],'dist':cat_vizer['Dist'],'Bmag':cat_vizer['Bmag'],'HyperLEDA':cat_vizer['HyperLEDA']})
            # ADD ME!!! #msk1=data[['dist']]<=distmax
            # ADD ME!!! #msk2=data[['dist']]>=distmin
            filter_good_ones =  (data.dist > 0) & \
                                  (data.dist !='NaN') & \
                                  (data.Bmag !='NaN') & \
                                  (data.Bmag !='null') & \
                                  (data.Bmag > 0)
            crossmatched_cat = data[filter_good_ones]
            outname = 'Xmatched_list_90_nside.csv'
            fullname = os.path.join(outdir, outname)
            crossmatched_cat.to_csv(fullname, sep="\t", index = False, header=True)
        if ORDERING == 'NUNIQ':
            start_moc_evt = time.time()
            uniq=data['UNIQ']
            probdensity=data['PROBDENSITY']
            level, ipix = ah.uniq_to_level_ipix(uniq)
            area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)
            prob = probdensity * area
            cumul_to = np.linspace(0.5, 0.9, 5)[::-1]
            colors = ['blue', 'green', 'yellow', 'orange', 'red']
            moxses = [MOC.from_valued_healpix_cells(uniq, prob, cumul_to=c) for c in cumul_to]
            df_gw['MOC']=moxses[0]
            print('moc_evt done and '+'Total Run-time is: '+str(time.time() - start_moc_evt))
            start_stmoc_evt = time.time()
            ###############
            #STMOC
            ###############
            outdir_stmoc = './GW_Alert/STMOC'
            if not os.path.exists(outdir_stmoc):
                os.mkdir(outdir_stmoc)
            string_list_gw=[]
            for i in range(0,len(df_gw)):
                string=str(df_gw['DATE_OBS'].values[i])
                string_list_gw.append(string)
            times_gw = Time(string_list_gw, format='isot', scale='utc')
            dt_iso_gw = TimeDelta(30*60, format='sec')
            t_GW_start=times_gw-dt_iso_gw
            t_GW_end=times_gw+dt_iso_gw
            df_gw['STMOC_GW_start']=t_GW_start
            df_gw['STMOC_GW_stop']=t_GW_end

            output_path='./GW_Alert/STMOC/AA_df_gw_STMOC_trial.csv'
            df_gw.to_csv(output_path, mode='a', header=not os.path.exists(output_path))

            stmoc_gw = STMOC.from_spatial_coverages(t_GW_start, t_GW_end, df_gw['MOC'])
            print("Time of the first observation: ", stmoc_gw.min_time.iso)
            print("Time of the last observation: ", stmoc_gw.max_time.iso)

            output_path='./GW_Alert/STMOC/'+str(graceid)+'_&_Revision_'+str(rev)+'_STMOC.fits'
            stmoc_gw.write(output_path,format='fits', overwrite=True)

            #Stacked
            #stacked_STMOC=STMOC.from_fits('./GW_Alert/STMOC/stacked_STMOC.fits')
            #restacked_STMOC=stmoc_gw.union(stacked_STMOC)
            #output_path='./GW_Alert/STMOC/stacked_STMOC.fits'
            #restacked_STMOC.write(output_path,format='fits', overwrite=True)

            if os.path.isfile('./GW_Alert/STMOC/AA_stacked_STMOC.fits'):
                print ("File does exist")
                stacked_STMOC=STMOC.from_fits('./GW_Alert/STMOC/AA_stacked_STMOC.fits')
                restacked_STMOC=stmoc_gw.union(stacked_STMOC)
                output_path='./GW_Alert/STMOC/AA_stacked_STMOC.fits'
                restacked_STMOC.write(output_path,format='fits', overwrite=True)
            else:
                print ("File does not exist")
                output_path='./GW_Alert/STMOC/AA_stacked_STMOC.fits'
                stmoc_gw.write(output_path,format='fits')
            print('stmoc_evt done and '+'Total Run-time is: '+str(time.time() - start_stmoc_evt))
            ###############
            #Figure NUNIQ
            ###############
            fig = plt.figure(111, figsize=(15, 10))
            with WCS(fig, 
                    fov=330 * u.deg,
                    center=SkyCoord(0,0, unit='deg', frame='galactic'),
                    coordsys='galactic',
                    rotation=Angle(0, u.degree),
                    projection="AIT") as wcs:
                ax = fig.add_subplot(1, 1, 1, projection=wcs)
                moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
                moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
                for (moc, c, col) in zip(moxses, cumul_to, colors):
                    moc.fill(ax=ax, wcs=wcs, alpha=0.5, linewidth=0, fill=True, color=col, label='confidence probability ' + str(round(c*100)) + '%')
                    moc.border(ax=ax, wcs=wcs, alpha=0.5, color=col)
                ax.legend()
            plt.xlabel('RA')
            plt.ylabel('DEC')
            plt.title('Bayestar')
            plt.grid(color="black", linestyle="dotted")
            outname = 'FOV_Galactic_nuniq.pdf'
            fullname = os.path.join(outdir, outname)    
            plt.savefig(fullname)
            print(len(moxses))
            print(cumul_to)
            moc_90_nuniq=moxses[0]
            cat_vizer = moc_90_nuniq.query_vizier_table('VII/281/glade2')
            data=pd.DataFrame({'RA': cat_vizer['_RAJ2000'], 'DEC':cat_vizer['_DEJ2000'],'dist':cat_vizer['Dist'],'Bmag':cat_vizer['Bmag'],'HyperLEDA':cat_vizer['HyperLEDA']})
            # ADD ME!!! #msk1=data[['dist']]<=distmax
            # ADD ME!!! #msk2=data[['dist']]>=distmin
            filter_good_ones =  (data.dist > 0) & \
                                  (data.dist !='NaN') & \
                                  (data.Bmag !='NaN') & \
                                  (data.Bmag !='null') & \
                                  (data.Bmag > 0)
            crossmatched_cat = data[filter_good_ones]
            outname = 'Xmatched_list_90_nuniq.csv'
            fullname = os.path.join(outdir, outname)    
            crossmatched_cat.to_csv(fullname, sep="\t", index = False, header=True)
        return crossmatched_cat, hdul1, outdir

    def Xmatched_raw_to_3Dplot(crossmatched_cat, hdul1, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat=crossmatched_cat.sort_values(by='Bmag', ascending=False)
        crossmatched_cat_top3=crossmatched_cat.head(3)
        header = hdul1[1].header
        DISTMEAN = header['DISTMEAN']
        DISTSTD = header['DISTSTD']
        DISTMIN = DISTMEAN-DISTSTD
        DISTMAX = DISTMEAN+DISTSTD
        ###############
        #Figure
        ###############
        fig = plt.figure(figsize=(20,10))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.view_init(20,-45)
        z = crossmatched_cat.dist.values
        x = crossmatched_cat.RA.values
        y = crossmatched_cat.DEC.values
        min_dist=min(z)
        ax.scatter(x, y, z, color='black', alpha=.9)
        X, Y = np.meshgrid(np.arange(min(x), max(x)), np.arange(min(y), max(y)))
        Z = 0*X+DISTMIN
        ax.plot_surface(X, Y, Z, color='b', alpha=.3)
        Z = 0*X+DISTMAX
        ax.plot_surface(X, Y, Z, color='g', alpha=.3)
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
        outname = 'Glade2_distribution.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)'''
        ########################################################
        ########################################################
        filter_distance =   (crossmatched_cat.dist <=DISTMAX) & \
                            (crossmatched_cat.dist >=DISTMIN)
        data_filtered_distance = crossmatched_cat[filter_distance]
        ########################################################
        ########################################################
        ###############
        #Figure
        ###############
        fig = plt.figure(figsize=(20,10))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.view_init(20,-45)
        z = data_filtered_distance.dist.values
        x = data_filtered_distance.RA.values
        y = data_filtered_distance.DEC.values
        min_dist=min(z)
        ax.scatter(x, y, z, color='black', alpha=.9)
        X, Y = np.meshgrid(np.arange(min(x), max(x)), np.arange(min(y), max(y)))
        Z = 0*X+DISTMIN
        ax.plot_surface(X, Y, Z, color='b', alpha=.3)
        Z = 0*X+DISTMAX
        ax.plot_surface(X, Y, Z, color='g', alpha=.3)
        m=ax.plot(x, y, 'ro', markersize=.5, color='r', zdir='z', zs=min_dist)
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.set_zlabel('Distance')
        ##################################################################################
        ##################################################################################
        N_bins = 30
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        ax.view_init(20,-45)
        hist, xedges, yedges = np.histogram2d(data_filtered_distance.RA.values, data_filtered_distance.DEC.values, bins=N_bins)
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
        outname = 'Glade2_Dist_filtered_catalog.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        return outdir

