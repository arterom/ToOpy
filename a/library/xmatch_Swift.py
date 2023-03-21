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
    def do_Xmatch(t_res, catalog, ra, dec, err, night, TrigID):
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
        outdir = './Swift_BAT_Alert/Xmatch/TRes'+str(t_res)+str('hrs')+'_&_'+str(observing_night)+'_&_TrigID_'+str(TrigID)
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
        idx_inside = moc_Swift_BAT.contains(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg)
        sources_inside = cat_4FGL[idx_inside]
        moc_Cross_matched = MOC.from_lonlat(sources_inside['_RAJ2000'].T * u.deg, sources_inside['_DEJ2000'].T * u.deg, 7)
        

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
        #Pandas df
        ###############
        dict_xmatch = {'RA': sources_inside['_RAJ2000'],
                'DEC': sources_inside['_DEJ2000'],
                '_4FGL': sources_inside['_4FGL'],
                'Cl1': sources_inside['Cl1'],
                'VarInd': sources_inside['VarInd'],
                'Assoc1': sources_inside['Assoc1']}
        crossmatched_cat=pd.DataFrame.from_dict(dict_xmatch)



        ###############
        #STMOC
        ###############
        outdir_stmoc = './STMOC'
        if not os.path.exists(outdir_stmoc):
            os.mkdir(outdir_stmoc)

        string_list_swift_bat=[]
        for i in range(0,len(df_swift)):
            string=str(df_swift['DATE-OBS'].values[i])
            string_list_swift_bat.append(string)
        times_bat = Time(string_list_swift_bat, format='isot', scale='utc')
        dt_iso_bat = TimeDelta(100, format='sec')
        t_bat_start=times_bat-dt_iso_bat
        t_bat_end=times_bat+dt_iso_bat
        df_swift['STMOC_bat_start']=t_bat_start
        df_swift['STMOC_bat_stop']=t_bat_end

        output_path='./STMOC/AA_df_swift_STMOC_trial.csv'
        df_swift.to_csv(output_path, mode='a', header=not os.path.exists(output_path))

        stmoc_bat = STMOC.from_spatial_coverages(t_bat_start, t_bat_end, df_swift['MOC'])
        print("Time of the first observation: ", stmoc_bat.min_time.iso)
        print("Time of the last observation: ", stmoc_bat.max_time.iso)

        output_path='./STMOC/Swift_TrigID_'+str(TrigID)+'_STMOC.fits'
        stmoc_bat.write(output_path,format='fits', overwrite=True)

        #Stacked
        #stacked_STMOC=STMOC.from_fits('./Swift_BAT_Alert/stacked_STMOC.fits')
        #restacked_STMOC=stmoc_bat.union(stacked_STMOC)
        #output_path='./Swift_BAT_Alert/stacked_STMOC.fits'
        #restacked_STMOC.write(output_path,format='fits', overwrite=True)

        if os.path.isfile('./STMOC/AA_stacked_STMOC.fits'):
            print ("File does exist")
            stacked_STMOC=STMOC.from_fits('./STMOC/AA_stacked_STMOC.fits')
            restacked_STMOC=stmoc_bat.union(stacked_STMOC)
            output_path='./STMOC/AA_stacked_STMOC.fits'
            restacked_STMOC.write(output_path,format='fits', overwrite=True)
        else:
            print ("File does not exist")
            output_path='./STMOC/AA_stacked_STMOC.fits'
            stmoc_bat.write(output_path,format='fits')



        ###############
        #Offsets
        ###############
        best_fit_positon = skycoord_evt
        CMatch = SkyCoord(crossmatched_cat.RA.values*u.deg, crossmatched_cat.DEC.values*u.deg)
        offset=best_fit_positon.separation(CMatch)/u.deg

        crossmatched_cat["offset"] = [float(item.value) for item in offset]
        crossmatched_cat=crossmatched_cat.sort_values(by='offset', ascending=True)

        ###############
        #Offset closest
        ###############
        if len(crossmatched_cat)==0:
            err_plus=float(err)+3
            moc_evt_enhanced = MOC.from_cone(
            lon=skycoord_evt.ra,
            lat=skycoord_evt.dec,
            radius=Angle(err_plus, u.deg),
            max_depth=10
            )
            idx_inside_enhanced = moc_evt_enhanced.contains(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg)
            sources_inside_enhanced = cat_4FGL[idx_inside_enhanced]
            moc_Cross_matched_enhanced = MOC.from_lonlat(sources_inside_enhanced['_RAJ2000'].T * u.deg, sources_inside_enhanced['_DEJ2000'].T * u.deg, 7)

            dict_xmatch_enhanced = {'RA': sources_inside_enhanced['_RAJ2000'],
                    'DEC': sources_inside_enhanced['_DEJ2000'],
                    '_4FGL': sources_inside_enhanced['_4FGL'],
                    'Cl1': sources_inside_enhanced['Cl1'],
                    'VarInd': sources_inside_enhanced['VarInd'],
                    'Assoc1': sources_inside_enhanced['Assoc1']}
            crossmatched_cat=pd.DataFrame.from_dict(dict_xmatch_enhanced)

            best_fit_positon = skycoord_evt
            CMatch = SkyCoord(crossmatched_cat.RA.values*u.deg, crossmatched_cat.DEC.values*u.deg)
            offset=best_fit_positon.separation(CMatch)/u.deg

            crossmatched_cat["offset"] = [float(item.value) for item in offset]
            crossmatched_cat=crossmatched_cat.sort_values(by='offset', ascending=True)
            crossmatched_cat=crossmatched_cat.head(1)
        ###############
        #Redshift Xmatch
        ###############
        #filter_Assoc_nan=(crossmatched_cat['Assoc1']=='')
        #df_filtered_nan = crossmatched_cat[filter_Assoc_nan]
        #df_filtered_nan["Redshift (from NED)"] = np.nan
        #df_filtered_nan["Assoc1"] = np.nan

        lst = ['','CRATES J023819+153323', 'Sim 147', 'Rosette', 'Monoceros', 'CRATES J081705+195836', 'CRATES J112431+230745', 'CRATES J112916+370317', 'CRATES J113514+301001']
        df_filtered_nan=crossmatched_cat.query('Assoc1 in @lst')
        df_filtered_nan["Redshift (from NED)"] = np.nan
        df_filtered_nan["Assoc1"] = np.nan
        
        df_filtered=crossmatched_cat.query('Assoc1 not in @lst')
        
        #filter_Assoc=(crossmatched_cat['Assoc1']!='')
        #df_filtered = crossmatched_cat[filter_Assoc]
        redshift_list=[]
        for i in range(0,len(df_filtered)):
            result_table = Ned.query_object(df_filtered['Assoc1'].values[i])
            redshift=result_table['Redshift']#[0]
            #print(df_filtered['Assoc1'].values[i])
            redshift_list.append(redshift)
        x=np.asarray(redshift_list)
        xnan = np.ma.filled(x.astype(float), np.nan)  ## Does not work anymore....
        #df_filtered["Redshift (from NED)"] = [item for item in xnan]
        df_filtered["Redshift (from NED)"] = [item for item in x]

        crossmatched_cat=pd.concat([df_filtered, df_filtered_nan])
        
        print('################################################################################################')
        print('##################### crossmatched_cat #########################################################')
        print('################################################################################################')
        print(crossmatched_cat)
        print('################################################################################################')
        print('################################################################################################')
        outname = 'Xmatched_list.csv'
        fullname = os.path.join(outdir, outname)    
        crossmatched_cat.to_csv(fullname, sep="\t", index = False, header=True)
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
        moc_Cross_matched.fill(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red", label='Xmatched')
        moc_Cross_matched.border(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red")
        moc_Swift_BAT.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Track Event')
        moc_Swift_BAT.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
        moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
        moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 20})
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
            center=SkyCoord(best_fit_positon.ra,best_fit_positon.dec, unit='deg', frame='icrs'),
            coordsys='icrs',
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc_Cross_matched.fill(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red", label='Xmatched')
        moc_Cross_matched.border(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red")
        moc_Swift_BAT.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Track Event')
        moc_Swift_BAT.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
        moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
        moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 20})
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
        return crossmatched_cat, crossmatched_cat_glade2, outdir

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


    def Xmatched_to_obslist(observatory, crossmatched_cat_glade2, zenith, moon_sep, time_resolution, night, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat_glade2=crossmatched_cat_glade2
        print('########################')
        print('########################')
        print(crossmatched_cat_glade2)
        print('########################')
        print('########################')
        crossmatched_cat_glade2=crossmatched_cat_glade2.sort_values(by='Bmag', ascending=False)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_swift.merged_def2.doit(observatory, crossmatched_cat_glade2, zenith, moon_sep, time_resolution, night, outdir)
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
