######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
#basics
import healpy as hp
import pandas as pd
import math
import os
import yaml
#astropy
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
from astropy.io import fits
from astropy.utils.data import download_file
from astropy.coordinates import Angle, SkyCoord
######################################################################################
######################################################################################
# Functions
######################################################################################
######################################################################################
class merged_def():
    def get_FermiLAT_weeks(event, rank):
        print('Output:')
        print('######################################################################################')
        print('######################################################################################')
        print('Bunch of warnings:')
        filename = download_file(event, cache=True)
        skymap_event, header = hp.read_map(filename, h=True, verbose=False)
        print('######################################################################################')
        print('######################################################################################')
        header = dict(header)
        NSIDE = header['NSIDE']
        EVENTID = header['EVENTID']
        START = header['START']
        print('Input of alert date (ISO) via header: '+str(START))
        EVENTMJD = header['EVENTMJD']
        print('Input of alert date (MJD) via header: '+str(EVENTMJD))
        print('######################################################################################')
        print('######################################################################################')
        EVENTMJD = int(EVENTMJD) #need otherwise weeks fail! fixit!
        print('Integer Input of alert date (MJD) via header: '+str(EVENTMJD))
        print('######################################################################################')
        print('######################################################################################')
        out_directory_date = header['START'][:11]
        I3TYPE = header['I3TYPE']
        RA=header['RA']
        print('RA location of the event (J2000 epoch): '+str(RA))
        DEC=header['DEC']
        print('DEC location of the event (J2000 epoch): '+str(DEC))
        CIRC_ERR90=header['CIRC_ERR90']
        print('Error90 location uncertainty (radius, statistical plus systematic, 90-percent containment): '+str(CIRC_ERR90))


        outdir = './Cascade_Alert/fermipy/test'+str(out_directory_date)
        if not os.path.exists(outdir):
            os.mkdir(outdir)   

        weeks_FermiLAT= [i for i in range(1,903)]
        df = pd.DataFrame(weeks_FermiLAT,columns =['week'])
        ##################################################################################
        ##################################################################################
        df["Iso_BL"] ='2008-06-05T00:00:00.000'
        df["MJD_BL"] ='54622.0'
        df["MET_BL"] =234316801.000
        ##################################################################################
        ##################################################################################
        df["dt_start"] = TimeDelta((df['week'].values-1)*7*86400, format='sec')
        df["dt_stop"] = TimeDelta((df['week'].values)*7*86400, format='sec')

        df["dt_start_mjd"] = TimeDelta((df['week'].values-1)*7, format='jd')
        df["dt_stop_mjd"] = TimeDelta((df['week'].values)*7, format='jd')

        df["dt_start_met"] = ((df['week'].values-1)*604800.0)
        df["dt_stop_met"] = ((df['week'].values)*604800.0)   
        ##################################################################################
        ##################################################################################
        df["Iso_Start"] = [Time(item, format='isot', scale='utc') for item in df["Iso_BL"].values] + df["dt_start"].values
        df["Iso_Stop"] = [Time(item, format='isot', scale='utc') for item in df["Iso_BL"].values] + df["dt_stop"].values

        df["Iso_Start_MJD"] = [Time(item, format='jd', scale='utc') for item in df["MJD_BL"].values] + df["dt_start_mjd"].values
        df["Iso_Stop_MJD"] = [Time(item, format='jd', scale='utc') for item in df["MJD_BL"].values] + df["dt_stop_mjd"].values

        df["Iso_Start_MET"] = df["MET_BL"].values + df["dt_start_met"].values
        df["Iso_Stop_MET"] = df["MET_BL"].values + df["dt_stop_met"].values
        ##################################################################################
        ##################################################################################
        df["Iso_Stop_MJD2"] = [item.jd for item in df["Iso_Stop_MJD"].values]
        df_event=df[df['Iso_Stop_MJD2'].between(EVENTMJD-3, EVENTMJD+3)]
        print('######################################################################################')
        print('######################################################################################')
        df_event=df_event.drop(columns=['Iso_BL', 'MJD_BL', 'MET_BL', 'dt_start', 'dt_stop', 'dt_start_mjd', 'dt_stop_mjd','dt_start_met', 'dt_stop_met', 'Iso_Stop_MJD2'])
        print(df_event)
        print('######################################################################################')
        print('######################################################################################')

        week_test=df_event['week'].values[0]
        print(week_test)
        outname = 'test '+str(EVENTID)+'& Date_'+str(START)+'.csv'
        fullname = os.path.join(outdir, outname)    
        df_event.to_csv(fullname, sep="\t", index = False, header=True)
        return header, df_event, outdir

    def prep_config_file(header, df_event, fermitools_refdata_path):
        header = dict(header)
        EVENTID = header['EVENTID']
        EVENTMJD = header['EVENTMJD']
        out_directory_date = header['START'][:11]
        RA=header['RA']
        DEC=header['DEC']
        CIRC_ERR90=header['CIRC_ERR90']
        #Event
        for i in range(0,len(df_event)):
            with open('../../../library/config_file_fermipy/config.yaml') as f:
                config = yaml.load(f, Loader=yaml.FullLoader)

            config['data']['evfile'] = './weekly/photon/lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits'
            config['data']['scfile'] = './weekly/spacecraft/lat_spacecraft_weekly_w'+str(df_event['week'].values[i])+'_p310_v001.fits'
            config['selection']['ra'] = RA
            config['selection']['dec'] = DEC
            #•  Time is expressed in UTC Mission Elapsed Time (MET)
            #•  You can convert MET into MJD (=modified Julian Date: MJD = JD – 2400000.5)
            config['selection']['tmin'] = int(df_event['Iso_Start_MET'].values[i])+4
            config['selection']['tmax'] = int(df_event['Iso_Stop_MET'].values[i])+4
            config['gtlike']['edisp'] = True
            config['model']['galdiff'] = str(fermitools_refdata_path)+'/fermi/galdiffuse/gll_iem_v07.fits'
            config['model']['isodiff'] = str(fermitools_refdata_path)+'/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt'
            
            with open('weekly_config.yaml', 'w') as f:
                config = yaml.dump(config, stream=f,
                               default_flow_style=False, sort_keys=False)