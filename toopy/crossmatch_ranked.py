######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
import argparse
import time
import os
import subprocess
from shlex import quote
from subprocess import call

from library import xmatch_track
from library import xmatch_cascade
from library import fermipy_methods_cascade
from library import fermipy_methods_track

from library import xmatch_GBM
from library import stmoc_GBM

from library import xmatch_Swift
from library import stmoc_Swift

from library import xmatch_GW
from library import tiled_GW



start_argparse = time.time()
######################################################################################
######################################################################################
# argparse
######################################################################################
######################################################################################
parser = argparse.ArgumentParser(description="Reads in arguments for each of the alert types")
######################################################################################
######################################################################################
# General
######################################################################################
######################################################################################
parser.add_argument("-observatory", "--observatory",
                    type=str,
                    help="Observatory performing the follow-up",
                    required=True)
parser.add_argument("-zenith", "--zenith",
                    type=float,
                    help="zenith constraint for observation",
                    required=True)
parser.add_argument("-moon_separation", "--moon_separation",
                    type=float,
                    help="moon_separation constraint for observation",
                    required=True)
parser.add_argument("-time_res", "--time_res",
                    type=float,
                    help="time resolution for observability",
                    required=True)
parser.add_argument("-flavour", "--flavour",
                    type=str,
                    help="Flavour of the Alert",
                    required=True,
                    choices=['Cascade', 'Track', 'GBM', 'GW', 'BAT'])
parser.add_argument("-mode", "--mode",
                    type=str,
                    help="mode",
                    required=False)
######################################################################################
######################################################################################
# FermiGBM & GW & CASCADE
######################################################################################
######################################################################################
parser.add_argument("-url", "--url",
                    type=str,
                    help="input url for GW event",
                    required=False)
parser.add_argument("-vol", "--vol",
                    type=float,
                    help="percent vol for Xmatch",
                    required=False)
######################################################################################
######################################################################################
# FermiGBM-ONLY
######################################################################################
######################################################################################
parser.add_argument("-TransNum_TrigID", "--TransNum_TrigID",
                    type=str,
                    help="Event distinction",
                    required=False)
######################################################################################
######################################################################################
# TRACK-ONLY
######################################################################################
######################################################################################
parser.add_argument("-ra", "--ra",
                    type=str,
                    help="EVENT RA (J2000) for Track Alert",
                    required=False)
parser.add_argument("-dec", "--dec",
                    type=str,
                    help="EVENT DEC (J2000) for Track Alert",
                    required=False)
parser.add_argument("-error", "--error",
                    type=str,
                    help="EVENT Error (90% CL, stat only) for Track Alert",
                    required=False)
parser.add_argument("-obs_night", "--obs_night",
                    type=str,
                    help="Observing night for Track Alert",
                    required=False)
parser.add_argument("-lightcurve", "--lightcurve",
                    type=str,
                    help="Boolean for lightcurve",
                    required=False)
######################################################################################
######################################################################################
# TRACK & CASCADE
######################################################################################
######################################################################################
parser.add_argument("-catalog", "--catalog",
                    type=str,
                    help="input catalog for crossmatch",
                    required=True)
parser.add_argument("-fermitools_refdata_path", "--fermitools_refdata_path",
                    type=str,
                    help="fermitools_refdata_path",
                    required=False)
parser.add_argument("-ranking_method", "--ranking_method",
                    type=str,
                    help="ranking methodology",
                    required=False,
                    default=False,
                    choices=['Xmatch', 'tiled_GW', 'STMOC', 'VarInd', 'FoV_prob', 'fermipy', 'etc'])
######################################################################################
######################################################################################
# TRACK & GW
######################################################################################
######################################################################################
parser.add_argument("-Rev", "--Rev",
                    type=str,
                    help="Revision of given alerts",
                    required=False)
######################################################################################
######################################################################################
# GW-ONLY
######################################################################################
######################################################################################
parser.add_argument("-GraceID", "--GraceID",
                    type=str,
                    help="Identifier in GraceDB",
                    required=False)
parser.add_argument("-instrument_FOV", "--instrument_FOV",
                    type=str,
                    help="Identifier in instrument_FOV",
                    required=False)
######################################################################################
######################################################################################
args = parser.parse_args()
# base directroy
######################################################################################
######################################################################################
base_dir=os.getcwd()
print(os.getcwd())
######################################################################################
######################################################################################
# Swift
######################################################################################
######################################################################################
if args.flavour == 'BAT':
    outdir = './Swift_BAT_Alert'
    if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
#######################################
#######################################
    if args.ranking_method == 'Xmatch':
        outdir = './Swift_BAT_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat_glade2, outdir=xmatch_Swift.merged_def.do_Xmatch(args.time_res, args.catalog, args.ra, args.dec, args.error, args.obs_night, args.TransNum_TrigID)
        outdir = xmatch_Swift.merged_def.Xmatched_raw_to_3Dplot(crossmatched_cat_glade2, outdir)
        xmatch_Swift.merged_def.Xmatched_to_obslist(args.observatory, crossmatched_cat_glade2, args.zenith, args.moon_separation, args.time_res, args.obs_night, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))

    if args.ranking_method == 'STMOC':
        outdir = './Swift_BAT_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat_glade2, outdir=stmoc_Swift.merged_def.do_Xmatch(args.time_res, args.catalog, args.ra, args.dec, args.error, args.obs_night, args.TransNum_TrigID)
        outdir = stmoc_Swift.merged_def.Xmatched_raw_to_3Dplot(crossmatched_cat_glade2, outdir)
        #xmatch_Swift.merged_def.Xmatched_to_obslist(args.observatory, crossmatched_cat_glade2, args.zenith, args.moon_separation, args.time_res, args.obs_night, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
        os.chdir(str(base_dir))
        argument_list=['test_arg']
        separator = " "
        subprocess.check_call('chmod u+r+x run_STMOC_Interceptor.sh', shell=True)
        subprocess.check_call("./run_STMOC_Interceptor.sh %s" % separator.join(argument_list), shell=True) 
######################################################################################
######################################################################################
# FermiGBM
######################################################################################
######################################################################################
if args.flavour == 'GBM':
    outdir = './GBM_Alert'
    if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
#######################################
#######################################
    if args.ranking_method == 'Xmatch':
        outdir = './GBM_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_glade2_cat, hdul1, outdir=xmatch_GBM.merged_def.do_Xmatch(args.url, args.TransNum_TrigID, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        outdir = xmatch_GBM.merged_def.Xmatched_raw_to_3Dplot_glade2(crossmatched_glade2_cat, hdul1, outdir)
        print('Done with xmatchGLade2 '+'Total Run-time is: '+str(time.time() - start_argparse))
        xmatch_GBM.merged_def.Xmatched_top10_BMag_to_obslist_glade2(args.url, args.observatory, crossmatched_glade2_cat, args.zenith, args.moon_separation, hdul1, args.time_res, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))

    if args.ranking_method == 'STMOC':
        outdir = './GBM_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_glade2_cat, hdul1, outdir=stmoc_GBM.merged_def.do_Xmatch(args.url, args.TransNum_TrigID, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        outdir = stmoc_GBM.merged_def.Xmatched_raw_to_3Dplot_glade2(crossmatched_glade2_cat, hdul1, outdir)
        print('Done with xmatchGLade2 '+'Total Run-time is: '+str(time.time() - start_argparse))
        #xmatch_GBM.merged_def.Xmatched_top10_BMag_to_obslist_glade2(args.url, args.observatory, crossmatched_glade2_cat, args.zenith, args.moon_separation, hdul1, args.time_res, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
        os.chdir(str(base_dir))
        argument_list=['test_arg']
        separator = " "
        subprocess.check_call('chmod u+r+x run_STMOC_Interceptor.sh', shell=True)
        subprocess.check_call("./run_STMOC_Interceptor.sh %s" % separator.join(argument_list), shell=True) 
######################################################################################
######################################################################################
# GW
######################################################################################
######################################################################################
if args.flavour == 'GW':
    outdir = './GW_Alert'
    if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
#######################################
#######################################
    if args.ranking_method == 'Xmatch':
        outdir = './GW_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, hdul1, outdir=xmatch_GW.merged_def.do_Xmatch(args.url, args.GraceID, args.Rev, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        outdir = xmatch_GW.merged_def.Xmatched_raw_to_3Dplot(crossmatched_cat, hdul1, outdir)
        xmatch_GW.merged_def.Xmatched_top10_BMag_to_obslist(args.url, args.observatory, crossmatched_cat, args.zenith, args.moon_separation, hdul1, args.time_res, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))

    if args.ranking_method == 'tiled_GW':
        outdir = './GW_Alert/tiled_GW'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, hdul1, outdir=tiled_GW.merged_def.do_Xmatch(args.url, args.GraceID, args.Rev, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog, args.mode, args.instrument_FOV)
        tiled_GW.merged_def.Xmatched_top10_tiling_to_obslist(args.url, args.observatory, crossmatched_cat, args.zenith, args.moon_separation, hdul1, args.time_res, args.mode, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
######################################################################################
######################################################################################
# TRACK
######################################################################################
######################################################################################
if args.flavour == 'Track':
    outdir = './Track_Alert'
    if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
    if args.ranking_method == 'VarInd':
        outdir = './Track_Alert/VarInd'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, outdir=xmatch_track.merged_def.do_Xmatch(args.time_res, args.catalog, args.ra, args.dec, args.error, args.obs_night, args.Rev)
        xmatch_track.merged_def.Xmatched_to_obslist(args.observatory, crossmatched_cat, args.zenith, args.moon_separation, args.obs_night, args.time_res, outdir)
        print('######################################################################################')
        print('######################################################################################')
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
        print('######################################################################################')
        print('######################################################################################')
##############################################################################################################################
##############################################################################################################################
    if args.ranking_method == 'fermipy':
        outdir = './Track_Alert/fermipy'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        df_event, outdir=fermipy_methods_track.merged_def.get_FermiLAT_weeks(args.ra, args.dec, args.error, args.obs_night, args.Rev)
        print(str(outdir))
        os.chdir(str(outdir))
##############################################################################################################################
##############################################################################################################################
        print('Checking whether FermiLat weekly data is available:')
        for i in range(0,len(df_event)):
            print('lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits')
            script_photon='wget -m -P . -nH --cut-dirs=4 -np -e robots=off https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits'
            script_spacecraft='wget -m -P . -nH --cut-dirs=4 -np -e robots=off https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/spacecraft/lat_spacecraft_weekly_w'+str(df_event['week'].values[i])+'_p310_v001.fits'
            PATH_photon = 'weekly/photon/lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits'
            if os.path.isfile(PATH_photon) and os.access(PATH_photon, os.R_OK):
                print("->Photon-file exists and is readable")
            else:
                print("!!!!!!!!!!->Either the spacecraft-file is missing or not readable<-!!!!!!!!!!")
                rc = call(script_photon, shell=True)
            PATH_spacecraft = 'weekly/spacecraft/lat_spacecraft_weekly_w'+str(df_event['week'].values[i])+'_p310_v001.fits'
            if os.path.isfile(PATH_spacecraft) and os.access(PATH_spacecraft, os.R_OK):
                print("->spacecraft-file exists and is readable")
            else:
                print("!!!!!!!!!!->Either the spacecraft-file is missing or not readable<-!!!!!!!!!!")
                rc = call(script_spacecraft, shell=True)
##############################################################################################################################
##############################################################################################################################
        config=fermipy_methods_track.merged_def.prep_config_file(args.ra, args.dec, df_event, args.obs_night, args.fermitools_refdata_path)
        print('######################################################################################')
        print('######################################################################################')
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
        print('######################################################################################')
        print('######################################################################################')
        print('Now follows the fermipy_runner script from bash:')
        print(str(base_dir))
##############################################################################################################################
##############################################################################################################################
        os.chdir(str(base_dir))
        argument_list=[str(outdir), args.lightcurve]
        separator = " "
        subprocess.check_call('chmod u+r+x run_fermipy.sh', shell=True)
        subprocess.check_call("./run_fermipy.sh %s" % separator.join(argument_list), shell=True) 
######################################################################################
######################################################################################
# CASCADE
######################################################################################
######################################################################################
if args.flavour == 'Cascade':
    outdir = './Cascade_Alert'
    if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
###############################################################
###############################################################
    if args.ranking_method == 'Xmatch':
        outdir = './Cascade_Alert/Xmatch'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, header, EVENTID, START, outdir=xmatch_cascade.merged_def.do_Xmatch(args.url, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        xmatch_cascade.merged_def.Xmatched_raw_to_obslist(crossmatched_cat, args.zenith, args.moon_separation, header, args.time_res, EVENTID, START, outdir)
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
###############################################################
###############################################################
    if args.ranking_method == 'VarInd':
        outdir = './Cascade_Alert/VarInd'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, header, EVENTID, START, outdir=xmatch_cascade.merged_def.do_Xmatch(args.url, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        xmatch_cascade.merged_def.Xmatched_top10_VarInd_to_obslist(args.url, args.observatory, crossmatched_cat, args.zenith, args.moon_separation, header, args.time_res, EVENTID, START, outdir)
        print('VarInd done and '+'Total Run-time is: '+str(time.time() - start_argparse))
###############################################################
###############################################################
    if args.ranking_method == 'FoV_prob':
        outdir = './Cascade_Alert/FoV_prob'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        crossmatched_cat, header, EVENTID, START, outdir=xmatch_cascade.merged_def.do_Xmatch(args.url, args.vol, args.ranking_method, args.time_res, args.zenith, args.catalog)
        xmatch_cascade.merged_def.Xmatched_top10_FoV_to_obslist(args.url, args.observatory, crossmatched_cat, args.zenith, args.moon_separation, header, args.time_res, EVENTID, START, outdir)
        print('FoV_prob done and '+'Total Run-time is: '+str(time.time() - start_argparse))
###############################################################
###############################################################
    if args.ranking_method == 'fermipy':
        outdir = './Cascade_Alert/fermipy'
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        header, df_event, outdir=fermipy_methods_cascade.merged_def.get_FermiLAT_weeks(args.url, args.ranking_method)
        print(str(outdir))
        os.chdir(str(outdir))
##############################################################################################################################
##############################################################################################################################
        print('Checking whether FermiLat weekly data is available:')
        for i in range(0,len(df_event)):
            print('lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits')
            script_photon='wget -m -P . -nH --cut-dirs=4 -np -e robots=off https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits'
            script_spacecraft='wget -m -P . -nH --cut-dirs=4 -np -e robots=off https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/spacecraft/lat_spacecraft_weekly_w'+str(df_event['week'].values[i])+'_p310_v001.fits'
            PATH_photon = 'weekly/photon/lat_photon_weekly_w'+str(df_event['week'].values[i])+'_p305_v001.fits'
            if os.path.isfile(PATH_photon) and os.access(PATH_photon, os.R_OK):
                print("->Photon-file exists and is readable")
            else:
                print("!!!!!!!!!!->Either the spacecraft-file is missing or not readable<-!!!!!!!!!!")
                rc = call(script_photon, shell=True)
            PATH_spacecraft = 'weekly/spacecraft/lat_spacecraft_weekly_w'+str(df_event['week'].values[i])+'_p310_v001.fits'
            if os.path.isfile(PATH_spacecraft) and os.access(PATH_spacecraft, os.R_OK):
                print("->spacecraft-file exists and is readable")
            else:
                print("!!!!!!!!!!->Either the spacecraft-file is missing or not readable<-!!!!!!!!!!")
                rc = call(script_spacecraft, shell=True)
##############################################################################################################################
##############################################################################################################################
        config=fermipy_methods_cascade.merged_def.prep_config_file(header, df_event, args.fermitools_refdata_path)
        print('######################################################################################')
        print('######################################################################################')
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
        print('######################################################################################')
        print('######################################################################################')
        print('Now follows the fermipy_runner script: ')
##############################################################################################################################
##############################################################################################################################
        os.chdir(str(base_dir))
        #subprocess.check_call('chmod u+r+x run_fermipy.sh', shell=True)
        #subprocess.check_call("./run_fermipy.sh %s" % quote(str(outdir)), shell=True)
        argument_list=[str(outdir), args.lightcurve]
        separator = " "
        subprocess.check_call('chmod u+r+x run_fermipy.sh', shell=True)
        subprocess.check_call("./run_fermipy.sh %s" % separator.join(argument_list), shell=True)
###############################################################
###############################################################
    if args.ranking_method == 'etc':
        print('-->!!Still need to work on that!!<--')
        print('Done and '+'Total Run-time is: '+str(time.time() - start_argparse))
