######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
import time
import gcn
import lxml.etree   
import subprocess
from urllib.request import urlopen
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.io import fits
from helper_gcnlistener.helper_functions import url_tester
from helper_gcnlistener.helper_functions import track_REV_tester
from helper_gcnlistener.helper_functions import find_between_tags
######################################################################################
######################################################################################
# Define your custom handler here.
######################################################################################
######################################################################################

# General Constraints
#############################################
#############################################
observatory='"Roque de los Muchachos"' 
max_zenith='50'
moon_separation='30'
too_span='daily' # scheduling time span to be performed for the trigger, i.e. daily, weekly or monthly.
time_resolution='1' # time in hours per observing slot

# Preferences for GW Alerts
#############################################
#############################################
instrument_FOV='15.5'
ranking='tiled_GW' # tiled_GW or Xmatch
mode='diagnostic' # diagnostic or performance

# Preferences for IceCube Alerts
#############################################
#############################################
fermitools_refdata_path='/opt/anaconda3/envs/toopy/share/fermitools/refdata'
lightcurve='no'


######################################################################################
#***To get all possible sites use astropy-module:
#from astropy.coordinates import EarthLocation
#astropy.coordinates.EarthLocation.get_site_names()
######################################################################################
######################################################################################
# Include notice_types
######################################################################################
######################################################################################
@gcn.include_notice_types(
    gcn.notice_types.SWIFT_BAT_GRB_POS_ACK,
    gcn.notice_types.FERMI_GBM_FLT_POS,
    gcn.notice_types.FERMI_GBM_GND_POS,
    gcn.notice_types.FERMI_GBM_FIN_POS,
    gcn.notice_types.FERMI_GBM_SUBTHRESH,
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.ICECUBE_ASTROTRACK_GOLD,
    gcn.notice_types.ICECUBE_ASTROTRACK_BRONZE,
    gcn.notice_types.ICECUBE_CASCADE)
######################################################################################
######################################################################################
# Handler
######################################################################################
######################################################################################
def handler(payload, root):
    parameter = {element.attrib['name']:
             element.attrib['value']
             for element in root.iterfind('./What/Param')}
    ###################################################################################### 
    #SWIFT_BAT_GRB_ALERT
    ######################################################################################
    ######################################################################################
    if '61' in parameter['Packet_Type']:
        print('we are looking at a SWIFT_BAT_GRB_POS_ACK Alert')
        print(root.attrib['ivorn'])
        TrigID=parameter['TrigID']

        outname='TrigID_'+str(TrigID)
        text_file = open('recieved_GCNs/SwiftBAT/SWIFT_BAT_GRB_POS_ACK_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()

        id=str(TrigID)
        url="https://gcn.gsfc.nasa.gov/other/%s.swift" % id
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        swift_url=url_tester(url, 30, 3)
        print('######################################################################################')
        print('######################################################################################')
        data = urlopen(swift_url).read()
        data=data.decode()
        text_file_webscrap = open('recieved_GCNs/SwiftBAT/SWIFT_BAT_GRB_POS_ACK_%s_webscrap.txt' % outname, 'w')
        text_file_webscrap.write(data)
        text_file_webscrap.close()
        print('#############################################')
        print('#############################################')
        print('Web-scrapping')
        print(data)
        for i, line in enumerate(data.split('\n')):
            if "GRB_RA" in line:
                RA=line
            if "GRB_DEC" in line:
                DEC=line
        RA_f=find_between_tags(RA,'{','}') 
        DEC_f=find_between_tags(DEC,'{','}')
        skycoord_evt = SkyCoord(RA_f,DEC_f, unit='deg', frame='icrs')
        EVT_ERROR90=3
        #ALERT TIME
        Event_TJD=parameter['Burst_TJD']
        Event_MJD=int(Event_TJD)+40000  # TJD = MJD - 40000
        Event_SOD=parameter['Burst_SOD']

        start_time_mjd = Time(Event_MJD, format='mjd', scale='utc')
        dt_sod = TimeDelta(Event_SOD, format='sec')
        start_time_mjd_actual = start_time_mjd + dt_sod
    ######################################################################################
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skycoord_evt.ra), str(skycoord_evt.dec), str(EVT_ERROR90), str(start_time_mjd_actual),  'Xmatch', str(TrigID)]
        separator = " "
        subprocess.check_call('chmod u+r+x ./method_scripts/SwiftBAT.sh', shell=True)
        subprocess.check_call("./method_scripts/SwiftBAT.sh %s" % separator.join(argument_list), shell=True)
    ######################################################################################
    ######################################################################################
    ###################################################################################### 
    #FERMI_GBM_FLT_POS
    ######################################################################################
    ######################################################################################
    if '111' in parameter['Packet_Type']:
        print('we are looking at a FERMI_GBM_FLT_POS Alert')
        print(root.attrib['ivorn'])
        TrigID=parameter['TrigID']

        outname='TrigID_'+str(TrigID)
        text_file = open('recieved_GCNs/FermiGBM/FLT-&GND-Pos/FERMI_GBM_FLT_POS_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()
        print('######################################################################################')
        print('######################################################################################')
        ######################################################################################
    ###################################################################################### 
    #FERMI_GBM_GND_POS
    ######################################################################################
    ######################################################################################
    if '112' in parameter['Packet_Type']:
        print('we are looking at a FERMI_GBM_GND_POS Alert')
        print(root.attrib['ivorn'])
        TrigID=parameter['TrigID']

        outname='TrigID_'+str(TrigID)
        text_file = open('recieved_GCNs/FermiGBM/FLT-&GND-Pos/FERMI_GBM_GND_POS_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()
        print('######################################################################################')
        print('######################################################################################')
    ######################################################################################
    ###################################################################################### 
    #FERMI_GBM_FIN_POS
    ######################################################################################
    ######################################################################################
    if '115' in parameter['Packet_Type']:
        print('we are looking at a FERMI_GBM_FIN_POS Alert')
        print(root.attrib['ivorn'])
        TrigID=parameter['TrigID']
        LocationMap_URL=parameter['LocationMap_URL']

        LocationMap_URL = str(LocationMap_URL)
        for r in (("locplot", "healpix"), ("png", "fit")):
            LocationMap_URL = LocationMap_URL.replace(*r)
        print(LocationMap_URL)

        outname='TrigID_'+str(TrigID)
        text_file = open('recieved_GCNs/FermiGBM/FERMI_GBM_FIN_POS_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()
        print('######################################################################################')
        print('######################################################################################')
        #These location maps are not available immediately -- the time delays are shown:
        #TYPE                 PKT_NUM     DELAY [min]
        #==================   =======     ===========
        #GBM_FIN_POS           115          1-15 min
        #https://gcn.gsfc.nasa.gov/fermi.html
        print('Attempting to reach url')
        skymap_fits=url_tester(LocationMap_URL, 3, 5)
        print('######################################################################################')
        print('######################################################################################')
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits), 'Xmatch', str(TrigID)+'_FinalPosition', too_span]
        separator = " "
        subprocess.check_call("./method_scripts/FermiGBM.sh %s" % separator.join(argument_list), shell=True)
    ######################################################################################
    ######################################################################################
    #FERMI_GBM_SUBTHRESH
    ######################################################################################
    ######################################################################################
    if '131' in parameter['Packet_Type']:
        print('we are looking at a FERMI_GBM_SUBTHRESH Alert')
        print(root.attrib['ivorn'])
        Trans_Num=parameter['Trans_Num']
        HealPix_URL=parameter['HealPix_URL']

        outname='Trans_Num_'+str(Trans_Num)
        text_file = open('recieved_GCNs/FermiGBM/FERMI_GBM_SUBTHRESH_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        skymap_fits=url_tester(HealPix_URL, 20, 3)
        print('######################################################################################')
        print('######################################################################################')
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits), 'Xmatch', str(Trans_Num)+'_SubTreshold', too_span]
        separator = " "
        subprocess.check_call("./method_scripts/FermiGBM.sh %s" % separator.join(argument_list), shell=True)
    ######################################################################################
    ######################################################################################
    #LVC_PRELIMINARY
    ######################################################################################
    ######################################################################################
    if '150' in parameter['Packet_Type']:
        print('we are looking at a LVC_PRELIMINARY Alert')
        GraceID=parameter['GraceID']
        Pkt_Ser_Num=parameter['Pkt_Ser_Num']
        Rev=int(Pkt_Ser_Num)-1
        outname='GraceID_'+str(GraceID)+'_Revision_'+str(Rev)
        text_file = open('recieved_GCNs/GW/LVC_PRELIMINARY_%s.txt' % outname, 'w')
        params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
        for key, value in params.items():
            print(key, '=', value)
            text_file.write('{} = {} \n'.format(key, value))
        text_file.close()
        if 'skymap_fits' in params:
            skymap_fits=params['skymap_fits']
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        skymap_fits_url=url_tester(skymap_fits, 30, 1)
        print('######################################################################################')
        print('######################################################################################')
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits_url), ranking, str(GraceID), str(Rev), mode, instrument_FOV]
        separator = " "
        subprocess.check_call("./method_scripts/GW.sh %s" % separator.join(argument_list), shell=True)
    ######################################################################################
    ######################################################################################
    #ICECUBE_ASTROTRACK_GOLD
    ######################################################################################
    ######################################################################################
    if '173' in parameter['Packet_Type']:
        print('we are looking at a ICECUBE_ASTROTRACK_GOLD Alert')
        print(root.attrib['ivorn'])
        RUN_ID=parameter['run_id']
        EVENT_ID=parameter['event_id']
        Rev=parameter['Rev']

        outname='RUN_ID_'+str(RUN_ID)+'_EVENT_ID_'+str(EVENT_ID)+'_Revision_'+str(Rev)
        text_file = open('recieved_GCNs/IceCube/GCN_ICECUBE_ASTROTRACK_GOLD_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()

        id=str(RUN_ID)+'_'+str(EVENT_ID)
        url="https://gcn.gsfc.nasa.gov/notices_amon_g_b/%s.amon" % id
        print(url)
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        Track_url=track_REV_tester(url, Rev, 30, 3)
        print('######################################################################################')
        print('######################################################################################')

        data = urlopen(Track_url).read()
        data=data.decode()
        text_file_webscrap = open('recieved_GCNs/IceCube/GCN_ICECUBE_ASTROTRACK_GOLD_%s_webscrap.txt' % outname, 'w')
        text_file_webscrap.write(data)
        text_file_webscrap.close()
        print('#############################################')
        print('#############################################')
        print('Web-scrapping')
        print(data)
        for i, line in enumerate(data.split('\n')):
            if "SRC_RA" in line:
                RA=line
            if "SRC_DEC" in line:
                DEC=line
        RA_f=find_between_tags(RA,'{','}') 
        DEC_f=find_between_tags(DEC,'{','}')
        skycoord_evt = SkyCoord(RA_f,DEC_f, unit='deg', frame='icrs')
        EVT_ERROR90=parameter['src_error_90']
        #ALERT TIME
        Event_TJD=parameter['Event_TJD']
        Event_MJD=int(Event_TJD)+40000  # TJD = MJD - 40000
        Event_SOD=parameter['Event_SOD']

        start_time_mjd = Time(Event_MJD, format='mjd', scale='utc')
        dt_sod = TimeDelta(Event_SOD, format='sec')
        start_time_mjd_actual = start_time_mjd + dt_sod
    ######################################################################################
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skycoord_evt.ra), str(skycoord_evt.dec), str(EVT_ERROR90), str(start_time_mjd_actual), 'VarInd', fermitools_refdata_path, str(Rev), 'no']
        separator = " "
        subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True) 
    ######################################################################################
        if '0' in parameter['Rev']:
            argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skycoord_evt.ra), str(skycoord_evt.dec), str(EVT_ERROR90), str(Event_MJD), 'fermipy', fermitools_refdata_path, str(Rev), lightcurve]
            separator = " "
            subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True)         
        if '1' in parameter['Rev']:
            print('No fermipy for you!') 
    ######################################################################################
    ######################################################################################
    #ICECUBE_ASTROTRACK_BRONZE
    ######################################################################################
    ######################################################################################
    if '174' in parameter['Packet_Type']:
        print('we are looking at a ICECUBE_ASTROTRACK_BRONZE Alert')
        print(root.attrib['ivorn'])
        RUN_ID=parameter['run_id']
        EVENT_ID=parameter['event_id']
        Rev=parameter['Rev']

        outname='RUN_ID_'+str(RUN_ID)+'_EVENT_ID_'+str(EVENT_ID)+'_Revision_'+str(Rev)
        text_file = open('recieved_GCNs/IceCube/GCN_ICECUBE_ASTROTRACK_BRONZE_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()

        id=str(RUN_ID)+'_'+str(EVENT_ID)
        url="https://gcn.gsfc.nasa.gov/notices_amon_g_b/%s.amon" % id
        print(url)
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        Track_url=track_REV_tester(url, Rev, 30, 3)
        print('######################################################################################')
        print('######################################################################################')
        data = urlopen(Track_url).read()
        data=data.decode()
        text_file_webscrap = open('recieved_GCNs/IceCube/GCN_ICECUBE_ASTROTRACK_BRONZE_%s_webscrap.txt' % outname, 'w')
        text_file_webscrap.write(data)
        text_file_webscrap.close()
        print('#############################################')
        print('#############################################')
        print('Web-scrapping')
        print(data)
        for i, line in enumerate(data.split('\n')):
            if "SRC_RA" in line:
                RA=line
            if "SRC_DEC" in line:
                DEC=line
        RA_f=find_between_tags(RA,'{','}') 
        DEC_f=find_between_tags(DEC,'{','}')
        skycoord_evt = SkyCoord(RA_f,DEC_f, unit='deg', frame='icrs')
        EVT_ERROR90=parameter['src_error_90']
        #ALERT TIME
        Event_TJD=parameter['Event_TJD']
        Event_MJD=int(Event_TJD)+40000  # TJD = MJD - 40000
        Event_SOD=parameter['Event_SOD']

        start_time_mjd = Time(Event_MJD, format='mjd', scale='utc')
        dt_sod = TimeDelta(Event_SOD, format='sec')
        start_time_mjd_actual = start_time_mjd + dt_sod
    ######################################################################################
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skycoord_evt.ra), str(skycoord_evt.dec), str(EVT_ERROR90), str(start_time_mjd_actual), 'VarInd', fermitools_refdata_path, str(Rev), 'no']
        separator = " "
        subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True) 
    ######################################################################################
        if '0' in parameter['Rev']:
            argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skycoord_evt.ra), str(skycoord_evt.dec), str(EVT_ERROR90), str(Event_MJD), 'fermipy', fermitools_refdata_path, str(Rev), lightcurve]
            separator = " "
            subprocess.check_call("./method_scripts/IceCube_TRACK.sh %s" % separator.join(argument_list), shell=True)         
        if '1' in parameter['Rev']:
            print('No fermipy for you!') 
    ######################################################################################
    ######################################################################################
    #ICECUBE_CASCADE
    ######################################################################################
    ######################################################################################
    if '176' in parameter['Packet_Type']:
        print('we are looking at a ICECUBE_CASCADE Alert')
        print(root.attrib['ivorn'])
        event_name=parameter['event_name']

        outname='Event_Name: '+str(event_name)
        text_file = open('recieved_GCNs/IceCube/GCN_ICECUBE_CASCADE_%s.txt' % outname, 'w')
        for param in root.findall('./What/Param'):
            name = param.attrib['name']
            value = param.attrib['value']
            print('{} = {}'.format(name, value))
            text_file.write('{} = {} \n'.format(name, value))
        text_file.close()

        skymap_fits=parameter['skymap_fits']
        print('######################################################################################')
        print('######################################################################################')
        print('Attempting to reach url')
        skymap_fits_url=url_tester(skymap_fits, 20, 3)
        print('######################################################################################')
        print('######################################################################################')
    ######################################################################################
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits_url), 'VarInd', fermitools_refdata_path, 'no', too_span]
        separator = " "
        subprocess.check_call("./method_scripts/IceCube_CASCADE.sh %s" % separator.join(argument_list), shell=True)
    ######################################################################################      
        argument_list=[observatory, max_zenith, moon_separation, time_resolution, str(skymap_fits_url), 'fermipy', fermitools_refdata_path, lightcurve, too_span]
        separator = " "
        subprocess.check_call("./method_scripts/IceCube_CASCADE.sh %s" % separator.join(argument_list), shell=True)

# Listen for VOEvents until killed with Control-C.
gcn.listen(handler=handler)

# Templates for test alerts
#payload = open('./xml_test_alerts/IceCube_CASCADE/2021/CASCADE_Dez21.xml', 'rb').read()
#root = lxml.etree.fromstring(payload)
#handler(payload, root)
