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
from astropy.table import vstack, Table
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
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
from library.helper import observability_gw_tiled
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(event, graceid, rev, vol_percent, rank, too_span, t_res, zenith, catalog, mode, instrument_FOV):
        ###############
        #Event
        ###############
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
        ###############
        #Rankings
        ###############
        if rank == 'tiled_GW':
            outdir = './GW_Alert/tiled_GW/'+str(too_span)+'_ToO_TRes'+str(t_res)+str('hrs')+'_&_'+str(out_directory_date)+'_&_'+str(graceid)+'_&_Revision_'+str(rev)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        ##################################
        #4FGL Catalog for Skymap
        ##################################
        vizier = VizierClass(
        row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
        cat_4FGL, = vizier.get_catalogs(catalog)
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        ##################################
        #Ra&Dec Grid for Skymap
        ##################################
        RA_list = [*range(0, 361, 2)]
        DEC_list = [*range(-90, 91, 2)]
        c = list(product(RA_list, DEC_list))
        RA = [x[0] for x in c]
        DEC = [x[1] for x in c]
        radec_list = pd.DataFrame(
                                 {'ra': RA,
                                  'dec': DEC
                                 })
        t2 = Table.from_pandas(radec_list)
        moc_radecgrid = MOC.from_lonlat(t2['ra'].T * u.deg, t2['dec'].T * u.deg, max_norder=7)
        print('######################################################################################')
        print('######################################################################################')
        print('############################Crab Reference##########################################')
        print('######################################################################################')
        FOV=float(instrument_FOV)
        d_crab = {'col0':'np.nan',
                        'RA': 83.6333,
                         'DEC': 22.0133,
                         'RA&DEC string': SkyCoord(83.6333*u.deg,22.0133*u.deg,frame='icrs').to_string('hmsdms'),
                         'Prob': 'np.nan',
                        'Tag':'Crab_check'}
        df_crab = pd.DataFrame(data=d_crab,index=[0])
        #spatial_coverage_Crab = MOC.from_cone(df_crab['RA'][0].T * u.deg,
        #                                      df_crab['DEC'][0].T * u.deg,
        #                                      FOV * u.deg, 9)
        #print(df_crab)
        ##################################
        #Most Probable Sky Location 
        #Following ->https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html#most-probable-sky-location
        ##################################    
        def Most_prob(event):
            hdul1 = fits.open(event)
            skymap = hdul1[1].data
            i = np.argmax(skymap['PROBDENSITY'])
            uniq = skymap[i]['UNIQ']
            rho_prob=skymap[i]['PROBDENSITY'] * (np.pi / 180)**2
            level, ipix = ah.uniq_to_level_ipix(uniq)
            nside = ah.level_to_nside(level)
            ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
            return ra.deg, dec.deg, rho_prob 
        first_pointing=Most_prob(event)
        d = {'RA': [first_pointing[0]],
             'DEC': [first_pointing[1]],
             'Prob': [first_pointing[2]],
             'Tag': ['1st_Pointing_Pix_centred']}
        df_c = pd.DataFrame(data=d)
        table_first = Table.from_pandas(df_c)
        positions_fc = SkyCoord(Most_prob(event)[0]* u.deg, Most_prob(event)[1]* u.deg)
        spatial_coverages_fc = MOC.from_cone(positions_fc.ra, positions_fc.dec, 3 * u.deg, 7)
        print('######################################################################################')
        print('######################################################################################')
        print('######################################################################################')
        print('######################################################################################')
        if ORDERING == 'NUNIQ':
            uniq=data['UNIQ']
            probdensity=data['PROBDENSITY']
            level, ipix = ah.uniq_to_level_ipix(uniq)
            area = ah.nside_to_pixel_area(ah.level_to_nside(level)).to_value(u.steradian)
            prob = probdensity * area
            cumul_to = np.linspace(0.5, 0.9, 5)[::-1]
            colors = ['blue', 'green', 'yellow', 'orange', 'red']
            moxses = [MOC.from_valued_healpix_cells(uniq, prob, cumul_to=c) for c in cumul_to]
            df_gw['MOC']=moxses[0]
            ##############################
            #Figure Skymap Galactic Coords
            ##############################
            if mode == 'diagnostic':
                fig = plt.figure(111, figsize=(15, 10))
                with WCS(fig, 
                        fov=360 * u.deg,
                        center=SkyCoord(0,0, unit='deg', frame='galactic'),
                        coordsys='galactic',
                        rotation=Angle(0, u.degree),
                        projection="AIT") as wcs:
                    ax = fig.add_subplot(1, 1, 1, projection=wcs)
                    moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
                    moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
                    #spatial_coverage_Crab.fill(ax=ax, wcs=wcs, alpha=0.1, color="red", fill=True, label='Pointing_on_Crab')
                    #spatial_coverage_Crab.border(ax=ax, wcs=wcs, alpha=0.5, fill=True)
                    #moc_radecgrid.fill(ax=ax, wcs=wcs, alpha=0.4, fill=True, color="red", label='RA & DEC Grid')
                    #moc_radecgrid.border(ax=ax, wcs=wcs, alpha=0.3, fill=True, color="black")
                    for (moc, c, col) in zip(moxses, cumul_to, colors):
                        moc.fill(ax=ax, wcs=wcs, alpha=0.5, linewidth=0, fill=True, color=col, label='confidence probability ' + str(round(c*100)) + '%')
                        moc.border(ax=ax, wcs=wcs, alpha=0.5, color=col)
                    ax.legend(prop={'size': 20}, loc='upper right')
                plt.xlabel('RA')
                plt.ylabel('DEC')
                plt.title('Bayestar')
                plt.grid(color="black", linestyle="dotted")
                outname = 'FOV_Galactic_nuniq.pdf'
                fullname = os.path.join(outdir, outname)    
                plt.savefig(fullname)
            print('######################################################################################')
            print('######################################################################################')
            print('######################################################################################')
            print('######################################################################################')
            moc_90p=moxses[0]
            idx_inside = moc_90p.contains_lonlat(t2['ra'].T * u.deg, t2['dec'].T * u.deg)
            sources_inside = t2[idx_inside]
            data_2=pd.DataFrame({
                                'RA': sources_inside['ra'], 
                               'DEC':sources_inside['dec'],
                               'RA&DEC string':'....................................',
                               'Prob':0})
            t3 = Table.from_pandas(data_2)
            moc_radecgrid_new = MOC.from_lonlat(t3['RA'].T * u.deg, t3['DEC'].T * u.deg, max_norder=7)
            ##################################
            #Most Probable Sky Location 
            #Following -> https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html#probability-density-at-a-known-position
            ##################################   
            def Prob_Galaxy(event,ra,dec):
                hdul1 = fits.open(event)
                skymap = hdul1[1].data
                level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
                nside = ah.level_to_nside(level)
                match_ipix = ah.lonlat_to_healpix(ra* u.deg, dec* u.deg, nside, order='nested')
                i = np.flatnonzero(ipix == match_ipix)[0]
                prob_density_4known_pos=skymap[i]['PROBDENSITY'] * (np.pi / 180)**2
                return ra, dec, prob_density_4known_pos
            for i in range(0,len(data_2)):
                ra, dec, prob_density_4known_pos=Prob_Galaxy(event, data_2.RA.values[i], data_2.DEC.values[i])
                data_2['Prob'][i]=prob_density_4known_pos
            data_2=data_2.sort_values(by=['Prob'],ascending=False)
            table_full = Table.from_pandas(data_2)
            df_complete=df_c.append(data_2)
            table_complete = Table.from_pandas(df_complete)
            print('######################################################################################')
            print('######################################################################################')
            print('############################Tiling Algorithm Start####################################')
            print('######################################################################################')
            print('######################################################################################')
            start_moc = moc_90p
            start_list = table_complete
            i=1
            pointing_list=[]
            while(start_moc.empty() == False):
                spatial_coverages_pointing = MOC.from_cone(start_list['RA'][0].T * u.deg,
                                             start_list['DEC'][0].T * u.deg,
                                             FOV * u.deg, 9)
                start_moc=start_moc.difference(spatial_coverages_pointing)
                pointing_list=vstack([pointing_list, start_list[0]])
                pointing_list['Tag']='Pointing_# '+str(i)
                ############################
                ############################
                # Figure for each Pointing
                ############################
                ############################
                if mode == 'diagnostic':
                    print('Working on pointing #'+str(i)+'...........')
                    fig = plt.figure(figsize=(15, 10))
                    with WCS(fig, 
                        fov=360 * u.deg,
                        center=SkyCoord(0,0, unit='deg', frame='galactic'),
                        coordsys='galactic',
                        rotation=Angle(0, u.degree),
                        projection="AIT") as wcs:
                        ax = fig.add_subplot(1, 1, 1, projection=wcs)
                    moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
                    moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
                    start_moc.fill(ax=ax, wcs=wcs, alpha=0.6, color='black', fill=True, label='To be covered')
                    start_moc.border(ax=ax, wcs=wcs, alpha=0.6, color='black', fill=True, label='')
                    spatial_coverages_pointing.fill(ax=ax, wcs=wcs, alpha=0.1, color="blue", fill=True, label='Pointing_# '+str(i))
                    spatial_coverages_pointing.border(ax=ax, wcs=wcs, alpha=0.5, fill=True)
                    #spatial_coverage_Crab.fill(ax=ax, wcs=wcs, alpha=0.1, color="red", fill=True, label='Pointing_on_Crab')
                    #spatial_coverage_Crab.border(ax=ax, wcs=wcs, alpha=0.5, fill=True)
                    moc_radecgrid_new.fill(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red", label='RA & DEC Grid')
                    moc_radecgrid_new.border(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="black")
                    ax.legend(prop={'size': 20}, loc='upper right')
                    plt.grid(color="black", linestyle="dotted")
                    outname = 'Pointing_# '+str(i)+'.pdf'
                    fullname = os.path.join(outdir, outname)    
                    plt.savefig(fullname)
                ####################
                ####################    
                idx_inside = start_moc.contains_lonlat(start_list['RA'].T * u.deg,
                                                start_list['DEC'].T * u.deg)
                sources_inside = start_list[idx_inside]
                i += 1
                if len(sources_inside) == 0:
                    break
                start_moc=start_moc.difference(spatial_coverages_pointing)
                start_list=start_list[idx_inside]
            print('######################################################################################')
            print('######################################################################################')
            print('############################Tiling Algorithm Finished#################################')
            print('######################################################################################')
            print('######################################################################################')
            ##################################
            #Cleaning up dataframe 
            ##################################  
            for i in range(0,len(pointing_list)):
                coords=SkyCoord(pointing_list['RA'][i]*u.deg,
                               pointing_list['DEC'][i]*u.deg,frame='icrs')
                pointing_list['Tag'][i]='Pointing_# '+str(i+1)
                pointing_list['RA&DEC string'][i]=coords.to_string('hmsdms')
            print(pointing_list)
            pointing_list_pre = pointing_list.to_pandas()

            crossmatched_cat=pd.concat([pointing_list_pre])#,df_crab])
            outname = 'tiling_list_90.csv'
            fullname = os.path.join(outdir, outname)    
            crossmatched_cat.to_csv(fullname, sep="\t", index = False, header=True)
        return crossmatched_cat, hdul1, outdir
    def Xmatched_top10_tiling_to_obslist(event, observatory, crossmatched_cat, zenith, moon_sep, hdul1, too_span, time_resolution, mode, outdir):
        ###############
        #Observability
        ###############
        ax, airmass_l, timetoplot, altitude_l, zenith_fin_l, c_fin, time_grid=observability_gw_tiled.merged_def2.doit(observatory, crossmatched_cat, zenith, moon_sep, hdul1, too_span, time_resolution, mode, outdir)
        print('######################################################################################')
        print('######################################################################################')
        print('############################Observability#############################################')
        print('######################################################################################')
        print('######################################################################################')
        ###############
        #Pandas
        ###############
        listed_obs=[]
        for i in range(0, len(c_fin)):
            dict = {'RA': crossmatched_cat.RA.values[i],
            'DEC': crossmatched_cat.DEC.values[i],
            'RA&DEC string':crossmatched_cat['RA&DEC string'].values[i],
            'Prob': crossmatched_cat.Prob.values[i],
            'Tag': crossmatched_cat.Tag.values[i],
            'Observable?': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Observatory': observatory,
            'Timeslot':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Airmass': [item for item in airmass_l[i]],
            'Altitude': [item/u.deg for item in altitude_l[i]],
            'Zenith': [item/u.deg for item in zenith_fin_l[i]]}  
            observability_df=pd.DataFrame(dict, index = [item for item in c_fin[i]])
            listed_obs.append(observability_df)
        listed_obs=pd.concat(listed_obs)
        outname = 'Observability_@ '+str(observatory)+'_AllSlots.csv'
        fullname = os.path.join(outdir, outname)    
        listed_obs.to_csv(fullname, sep="\t", index = False, header=True)
        print(listed_obs)
        print('This is zenith:'+str(zenith))
        print('This is time_resolution:'+str(time_resolution))
        fin_df=listed_obs.loc[True]
        outname = 'Observability_@ '+str(observatory)+'_TrueFlag.csv'
        fullname = os.path.join(outdir, outname)    
        fin_df.to_csv(fullname, sep="\t", index = False, header=True)
