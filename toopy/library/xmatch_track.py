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
from mocpy import WCS
#vizier
from astroquery.vizier import VizierClass
from astroquery.ned import Ned
#helper
from library.helper import observability_track
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(t_res, catalog, ra, dec, err, night, rev):
        ###############
        #Event
        ###############
        observing_night=Time(night, format='mjd', scale='utc').isot 
        print('Manual input of alert date: '+str(observing_night))
        outdir = './Track_Alert/VarInd/TRes'+str(t_res)+str('hrs')+'_&_'+str(observing_night)+'_&_Revision_'+str(rev)
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
        moc_Gold_n_Bronze = MOC.from_cone(
        lon=skycoord_evt.ra,
        lat=skycoord_evt.dec,
        radius=Angle(err, u.deg),
        max_depth=10
        )
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        idx_inside = moc_Gold_n_Bronze.contains(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg)
        sources_inside = cat_4FGL[idx_inside]
        moc_Cross_matched = MOC.from_lonlat(sources_inside['_RAJ2000'].T * u.deg, sources_inside['_DEJ2000'].T * u.deg, 7)
        
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
        moc_Gold_n_Bronze.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Track Event')
        moc_Gold_n_Bronze.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
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
        moc_Gold_n_Bronze.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Track Event')
        moc_Gold_n_Bronze.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
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
        return crossmatched_cat, outdir

    def Xmatched_to_obslist(observatory, crossmatched_cat, zenith, moon_sep, time_resolution, night, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat=crossmatched_cat
        crossmatched_cat=crossmatched_cat.sort_values(by='VarInd', ascending=False)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_track.merged_def2.doit(observatory, crossmatched_cat, zenith, moon_sep, time_resolution, night, outdir)
        ###############
        #Pandas
        ###############
        listed_obs=[]
        for i in range(0, len(c_fin)):
            dict = {'RA': crossmatched_cat.RA.values[i],
            'DEC': crossmatched_cat.DEC.values[i],
            '_4FGL': crossmatched_cat._4FGL.values[i],
            'Cl1': crossmatched_cat.Cl1.values[i],
            'VarInd': crossmatched_cat.VarInd.values[i],
            'Assoc1': crossmatched_cat.Assoc1.values[i],
            'Redshift': crossmatched_cat['Redshift (from NED)'].values[i],
            'Observable': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Observatory': observatory,
            'Timeslot':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Airmass': airmass,
            'Altitude': altitude,
            'Zenith': zenith,
            'offset': crossmatched_cat.offset.values[i]}
            observability_df=pd.DataFrame(dict, index = [item for item in c_fin[i]])
            listed_obs.append(observability_df)
        listed_obs=pd.concat(listed_obs)
        fin_df=listed_obs.loc[True]
        outname = 'Observability_@'+str(observatory)+'.csv'
        fullname = os.path.join(outdir, outname)    
        fin_df.to_csv(fullname, sep="\t", index = False, header=True)
