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
import os
#astropy
import astropy.units as u
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.coordinates import Angle, SkyCoord
#mocpy
from mocpy import MOC
from mocpy import WCS
#vizier
from astroquery.vizier import VizierClass
from astroquery.ned import Ned
#helper
from library.helper import observability_cascade
from library.helper import helper_NED_query
######################################################################################
######################################################################################
#Deal with "urllib.error.URLError: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: certificate has expired (_ssl.c:1091)>
######################################################################################
######################################################################################
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def do_Xmatch(event, vol_percent, rank, t_res, zenith, catalog):
        ###############
        #Event
        ###############
        filename = download_file(event, cache=True, allow_insecure=True)
        skymap_event, header = hp.read_map(filename, h=True, verbose=False)
        quantile = vol_percent # to get 50% contour
        header = dict(header)
        NSIDE = header['NSIDE']
        EVENTID = header['EVENTID']
        print(header)
        START = header['START']
        RA_header = header['RA']
        DEC_header = header['DEC']
        CIRC_ERR90 = header['CIRC_ERR90']
        print('This is START:'+str(START))
        EVENTMJD = header['EVENTMJD']
        out_directory_date = header['START'][:11]
        I3TYPE = header['I3TYPE']

        if rank == 'Xmatch':
            outdir = './Cascade_Alert/Xmatch/TRes'+str(t_res)+str('hrs')+'_&_'+str(out_directory_date)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        if rank == 'VarInd':
            outdir = './Cascade_Alert/VarInd/TRes'+str(t_res)+str('hrs')+'_&_'+str(out_directory_date)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        if rank == 'FoV_prob':
            outdir = './Cascade_Alert/FoV_prob/TRes'+str(t_res)+str('hrs')+'_&_'+str(out_directory_date)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
        #outdir = './Ranking_Method='+str(rank)+'Time_res='+str(t_res)+str('hrs')+str(vol_percent*100)+str('%')+'Event_type='+str(I3TYPE)+str('_Date_')+str(out_directory_date)+str('_Volume_')
        #if not os.path.exists(outdir):
        #    os.mkdir(outdir)

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
        '''
        #MAGIC exposure
        exposure = hpt.exposure_pdf(64, a0=28.761944, zmax=zenith)
        quantile = 1 # to get full contour
        argsort = np.argsort(-exposure)
        cum_skymap_exp = np.cumsum(sorted(exposure,reverse=True))
        cont_ind = argsort[cum_skymap_exp < quantile]
        contour = np.array([1. if pix in cont_ind else 0. for pix in range(len(exposure))])
        max_pix = np.argmax(contour)
        wh_c=np.where(contour == 1)[0]
        wh2_c=np.argwhere(contour == 1)
        dec, ra = hp.pix2ang(nside=64, ipix=[max_pix])
        dec = math.pi/2. - dec
        theta, phi = hp.pix2ang(64, wh_c)
        ra_Mexp = np.rad2deg(phi)
        dec_Mexp = np.rad2deg(0.5 * np.pi - theta)
        skycoord_Mexp = SkyCoord(ra_Mexp, dec_Mexp, unit="deg", frame="icrs")
        '''
        ###############
        #4FGL
        ###############
        vizier = VizierClass(
        row_limit=-1, columns=[ '*', '_RAJ2000', '_DEJ2000'])
        cat_4FGL, = vizier.get_catalogs(catalog)
        ###############
        #MOC
        ###############
        # With MAGIC_exposure, still fails
        '''
        moc_Cascade_evt = MOC.from_skycoords(skycoord_Cascade_evt, max_norder=7)
        moc_Mexp = MOC.from_skycoords(skycoord_Mexp, max_norder=6)
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        #Cross-match
        #4FGL&MAGIC_exp
        idx_inside = moc_Mexp.contains(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg)
        sources_inside = cat_4FGL[idx_inside]
        #Uncertainty&4FGL
        idx2_inside = moc_Cascade_evt.contains(sources_inside['_RAJ2000'].T * u.deg, sources_inside['_DEJ2000'].T * u.deg)
        sources2_inside = sources_inside[idx2_inside]
        moc_Cross_matched = MOC.from_lonlat(sources2_inside['_RAJ2000'].T * u.deg, sources2_inside['_DEJ2000'].T * u.deg, 7)
        #PANDAS

        dict_xmatch = {'RA': sources2_inside['_RAJ2000'],
                'DEC': sources2_inside['_DEJ2000'],
                '_4FGL': sources2_inside['_4FGL'],
                'Cl1': sources2_inside['Cl1'],
                'VarInd': sources2_inside['VarInd'],
                'Assoc1': sources2_inside['Assoc1']}
        crossmatched_cat=pd.DataFrame.from_dict(dict_xmatch)
        '''
        # Without MAGIC_exposure, still fails
        moc_Cascade_evt = MOC.from_skycoords(skycoord_Cascade_evt, max_norder=7)
        #moc_Mexp = MOC.from_skycoords(skycoord_Mexp, max_norder=6)
        moc_4FGL = MOC.from_lonlat(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg, max_norder=7)
        ###############
        #XMatch
        ###############
        idx_inside = moc_Cascade_evt.contains(cat_4FGL['_RAJ2000'].T * u.deg, cat_4FGL['_DEJ2000'].T * u.deg)
        sources_inside = cat_4FGL[idx_inside]
        print('This is len of sources_inside: '+str(len(sources_inside)))
        #Uncertainty&4FGL
        moc_Cross_matched = MOC.from_lonlat(sources_inside['_RAJ2000'].T * u.deg, sources_inside['_DEJ2000'].T * u.deg, 7)
        ###############
        #Pandas
        ###############
        dict_xmatch = {'RA': sources_inside['_RAJ2000'],
                'DEC': sources_inside['_DEJ2000'],
                '_4FGL': sources_inside['_4FGL'],
                'Cl1': sources_inside['Cl1'],
                'VarInd': sources_inside['VarInd'],
                'Assoc1': sources_inside['Assoc1']}
        crossmatched_cat=pd.DataFrame.from_dict(dict_xmatch)
        '''
        for col, dtype in crossmatched_cat.dtypes.items():
            if dtype == np.object:  # Only process byte object columns.
                crossmatched_cat[col] = crossmatched_cat[col].apply(lambda x: x.decode("utf-8"))
        '''
        ###############
        #Offsets
        ###############
        skycoord_evt = SkyCoord(RA_header, DEC_header, unit='deg', frame='icrs')
        best_fit_positon = skycoord_evt
        CMatch = SkyCoord(crossmatched_cat.RA.values*u.deg, crossmatched_cat.DEC.values*u.deg)
        offset=best_fit_positon.separation(CMatch)/u.deg
        crossmatched_cat["offset"] = [float(item.value) for item in offset]
        #crossmatched_cat=crossmatched_cat.sort_values(by='offset', ascending=True)
        ###############
        #Redshift
        ###############
        lst = ['','CRATES J023819+153323', 'Sim 147', 'Rosette', 'Monoceros', 'CRATES J081705+195836', 'CRATES J112431+230745', 'CRATES J112916+370317', 'CRATES J113514+301001', 'CRATES J060650+440144', 'LMC-FarWest', 'LMC SNR N 132D', 'LMC-30DorWest', 'LMC-North', 'LMC P3']
        df_filtered_nan=crossmatched_cat.query('Assoc1 in @lst')
        df_filtered_nan["Redshift (from NED)"] = np.nan
        df_filtered_nan["Assoc1"] = np.nan
        df_filtered=crossmatched_cat.query('Assoc1 not in @lst')
        #df_filtered["Redshift (from NED)"] = np.nan
        redshift_list=[]
        for i in range(0,len(df_filtered)):
            #print(df_filtered['Assoc1'].values[i])
            #result_table = Ned.query_object(df_filtered['Assoc1'].values[i])
            #redshift=result_table['Redshift']#[0]
            redshift=helper_NED_query.query_NED_object(df_filtered['Assoc1'].values[i])
            #print(redshift)
            redshift_list.append(redshift)
        #print(redshift_list)
        #x=np.asarray(redshift_list)
        #print(x)
        #xnan = np.ma.filled(x.astype(float), np.nan) ## Does not work anymore....
        #df_filtered["Redshift (from NED)"] = [item for item in xnan]
        df_filtered["Redshift (from NED)"] = [item for item in redshift_list]
        
        crossmatched_cat=pd.concat([df_filtered, df_filtered_nan])
        ###############
        #Offset-Sorting
        ###############        
        crossmatched_cat=crossmatched_cat.sort_values(by='offset', ascending=True)
        ###############
        #Print dataframe
        ############### 
        print('#############################################################################################################')
        print('################################## crossmatched_cat #########################################################')
        print('#############################################################################################################')
        print(crossmatched_cat)
        print('#############################################################################################################')
        print('#############################################################################################################')
        outname = 'Xmatched_list.csv'
        fullname = os.path.join(outdir, outname)    
        crossmatched_cat.to_csv(fullname, sep="\t", index = False, header=True)
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
            moc_Cross_matched.fill(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red", label='Xmatched')
            moc_Cross_matched.border(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red")
            moc_Cascade_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Cascade ('+str(vol_percent*100)+str('%')+' containment)')
            moc_Cascade_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
            #moc_Mexp.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='MAGIC Exposure 1 yr (WIP)')
            #moc_Mexp.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
            moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 10})
        plt.grid(color="black", linestyle="dotted")
        outname = 'FOV_Galactic.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        
        fig = plt.figure(figsize=(10, 10))
        with WCS(fig, 
            fov=50 * u.deg,
            center=SkyCoord(RA_header,DEC_header, unit='deg', frame='icrs'),
            coordsys='icrs',
            rotation=Angle(0, u.degree),
            projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc_Cross_matched.fill(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red", label='Xmatched')
        moc_Cross_matched.border(ax=ax, wcs=wcs, alpha=0.8, fill=True, color="red")
        moc_Cascade_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Cascade ('+str(vol_percent*100)+str('%')+' containment)')
        moc_Cascade_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
        #moc_Mexp.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='MAGIC Exposure 1 yr (WIP)')
        #moc_Mexp.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
        moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
        ax.legend(prop={'size': 20})
        plt.grid(color="black", linestyle="dotted")
        plt.xlabel('RA', size=30)
        plt.ylabel('DEC', size=30)
        plt.xticks(color='black', fontsize=22)
        ax.coords.grid(True, color='black')
        ax.coords[0].set_ticklabel(size="xx-large")
        ax.coords[1].set_ticklabel(size="xx-large")
        outname = 'FOV_ROI.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
        
        return crossmatched_cat, header, EVENTID, START,  outdir
    def Xmatched_raw_to_obslist(crossmatched_cat, zenith, moon_sep, header, time_resolution, EVENTID, START, outdir):
        ###############
        #Event
        ###############
        crossmatched_cat=crossmatched_cat
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_cascade.merged_def2.doit(crossmatched_cat, zenith, moon_sep, header, time_resolution, outdir)
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
            'Observable?': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Timeslot':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Airmass': airmass,
           'Altitude': altitude,
            'Zenith': zenith,
            'offset': crossmatched_cat.offset.values[i]}
            observability_df=pd.DataFrame(dict, index = [item for item in c_fin[i]])
            listed_obs.append(observability_df)
        listed_obs=pd.concat(listed_obs)
        fin_df=listed_obs.loc[True]
        outname = 'Observability_listed.csv'
        fullname = os.path.join(outdir, outname)    
        fin_df.to_csv(fullname, sep="\t", index = False, header=True)
    def Xmatched_top10_VarInd_to_obslist(event, observatory, crossmatched_cat, zenith, moon_sep, header,  time_resolution, EVENTID, START, outdir):
        ###############
        #Event
        ###############
        filename = download_file(event, cache=True)
        skymap_event, header = hp.read_map(filename, h=True, verbose=False)
        def circ_prob(ra, dec, radius):
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            radius = np.deg2rad(radius)
            xyz = hp.ang2vec(theta, phi)
            ipix_disc = hp.query_disc(128, xyz, radius)
            summed=skymap_event[ipix_disc].sum()
            return summed
        listofpercent_all_together=[]
        for i in range(0, len(crossmatched_cat)):
            summed=circ_prob(crossmatched_cat.RA.values[i], crossmatched_cat.DEC.values[i], 3.5)
            listofpercent_all_together.append(summed)
        crossmatched_cat['FoV_Prob'] = listofpercent_all_together
        crossmatched_cat=crossmatched_cat.sort_values(by='VarInd', ascending=False)
        crossmatched_cat=crossmatched_cat.head(10)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_cascade.merged_def2.doit(observatory, crossmatched_cat, zenith, moon_sep, header, time_resolution, outdir)
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
            'FoV_Prob': crossmatched_cat.FoV_Prob.values[i],
            'Assoc1': crossmatched_cat.Assoc1.values[i],
            'Redshift': crossmatched_cat['Redshift (from NED)'].values[i],
            'Observatory': observatory,
            'Observable?': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
            'Timeslot Start':[t.datetime.strftime("%H:%M") for t in time_grid],
            'Timeslot Duration [hr]':time_resolution,
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
        ###############
        #Figure
        ###############
        positions = SkyCoord(crossmatched_cat['RA'].values*u.deg, crossmatched_cat['DEC'].values*u.deg)
        spatial_coverages = [MOC.from_cone(pos.ra, pos.dec, 3.5 * u.deg, 7) for pos in positions]
        fig = plt.figure(figsize=(5,5))
        with WCS(fig, 
                fov=360 * u.deg,
                #center=SkyCoord(positions.ra[0],positions.dec[0], unit='deg', frame='galactic'),
                center=SkyCoord(0,0, unit='deg', frame='galactic'),
                coordsys='galactic',
                rotation=Angle(0, u.degree),
                projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
            '''
            moc_Cascade_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Cascade ('+str(vol_percent*100)+str('%')+' containment)')
            moc_Cascade_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
            moc_Mexp.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='MAGIC Exposure 1 yr')
            moc_Mexp.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
            moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            '''
            for i in range(0, len(spatial_coverages)):
                spatial_coverages[i].fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, label=''+str(crossmatched_cat['_4FGL'].values[i])+'')
                spatial_coverages[i].border(ax=ax, wcs=wcs, alpha=0.5, fill=True)
        ax.legend(prop={'size': 20})
        plt.grid(color="black", linestyle="dotted")
        outname = 'Observability_listed.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
    def Xmatched_top10_FoV_to_obslist(event, observatory, crossmatched_cat, zenith, moon_sep, header, time_resolution, EVENTID, START, outdir):
        ###############
        #Event
        ###############
        filename = download_file(event, cache=True)
        skymap_event, header = hp.read_map(filename, h=True, verbose=False)
        #skymap, header = hp.read_map(filename, h=True, verbose = True)
        def circ_prob(ra, dec, radius):
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            radius = np.deg2rad(radius)
            xyz = hp.ang2vec(theta, phi)
            ipix_disc = hp.query_disc(128, xyz, radius)
            summed=skymap_event[ipix_disc].sum()
            return summed

        listofpercent_all_together=[]
        for i in range(0, len(crossmatched_cat)):
            summed=circ_prob(crossmatched_cat.RA.values[i], crossmatched_cat.DEC.values[i], 3.5)
            listofpercent_all_together.append(summed)
        crossmatched_cat['FoV_Prob'] = listofpercent_all_together
        crossmatched_cat=crossmatched_cat.sort_values(by='FoV_Prob', ascending=False)
        crossmatched_cat=crossmatched_cat.head(10)
        ###############
        #Observability
        ###############
        ax, airmass, timetoplot, altitude, zenith, c_fin, time_grid=observability_cascade.merged_def2.doit(observatory, crossmatched_cat, zenith, moon_sep, header, time_resolution, outdir)
        ###############
        #Pandas
        ###############
        print([t.datetime.strftime("%H:%M") for t in time_grid])
        listed_obs=[]
        for i in range(0, len(c_fin)):
            dict = {'RA': crossmatched_cat.RA.values[i],
            'DEC': crossmatched_cat.DEC.values[i],
            '_4FGL': crossmatched_cat._4FGL.values[i],
            'Cl1': crossmatched_cat.Cl1.values[i],
            'VarInd': crossmatched_cat.VarInd.values[i],
            'FoV_Prob': crossmatched_cat.FoV_Prob.values[i],
            'Assoc1': crossmatched_cat.Assoc1.values[i],
            'Redshift': crossmatched_cat['Redshift (from NED)'].values[i],
            'Observatory': observatory,
            'Observable?': [item for item in c_fin[i]],
            'Observing Night': [t.datetime.strftime("%D") for t in time_grid],
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
        ###############
        #Figure
        ###############
        positions = SkyCoord(crossmatched_cat['RA'].values*u.deg, crossmatched_cat['DEC'].values*u.deg)
        spatial_coverages = [MOC.from_cone(pos.ra, pos.dec, 3.5 * u.deg, 7) for pos in positions]
        fig = plt.figure(figsize=(5,5))
        with WCS(fig, 
                fov=360 * u.deg,
                #center=SkyCoord(positions.ra[0],positions.dec[0], unit='deg', frame='galactic'),
                center=SkyCoord(0,0, unit='deg', frame='galactic'),
                coordsys='galactic',
                rotation=Angle(0, u.degree),
                projection="AIT") as wcs:
            ax = fig.add_subplot(1, 1, 1, projection=wcs)
            '''
            moc_Cascade_evt.fill(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="grey", label='Cascade ('+str(vol_percent*100)+str('%')+' containment)')
            moc_Cascade_evt.border(ax=ax, wcs=wcs, alpha=0.2, fill=True, color="black")
            moc_Mexp.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="green", label='MAGIC Exposure 1 yr')
            moc_Mexp.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            moc_4FGL.fill(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="blue", label='4FGL')
            moc_4FGL.border(ax=ax, wcs=wcs, alpha=0.1, fill=True, color="black")
            '''
            for i in range(0, len(spatial_coverages)):
                spatial_coverages[i].fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, label=''+str(crossmatched_cat['_4FGL'].values[i])+'')
                spatial_coverages[i].border(ax=ax, wcs=wcs, alpha=0.5, fill=True)
        ax.legend(prop={'size': 20})
        plt.grid(color="black", linestyle="dotted")
        outname = 'Observability_listed.pdf'
        fullname = os.path.join(outdir, outname)    
        plt.savefig(fullname)
