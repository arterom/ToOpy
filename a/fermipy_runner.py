######################################################################################
######################################################################################
# Imports
######################################################################################
######################################################################################
import argparse
import time
import yaml
import os
import pandas as pd
from fermipy.gtanalysis import GTAnalysis
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter
import matplotlib.pyplot as plt
from astropy.time import Time
plt.switch_backend('agg')
######################################################################################
######################################################################################
# Args
######################################################################################
######################################################################################
start_fermipy = time.time()
parser = argparse.ArgumentParser(description="etc")
parser.add_argument("-outpath", "--outpath",
                    type=str,
                    help="outpath",
                    required=True)
parser.add_argument("-lightcurve", "--lightcurve",
                    type=str,
                    help="lightcurve",
                    required=True)
args = parser.parse_args()
os.chdir(args.outpath)
######################################################################################
######################################################################################
# Fermipy Analysis
######################################################################################
######################################################################################
gta = GTAnalysis('weekly_config.yaml',logging={'verbosity' : 3})

with open('weekly_config.yaml') as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
tmin=config['selection']['tmin']
tmax=config['selection']['tmax']

gta.setup()
print('######################################################################################')
print('######################################################################################')
print('#################### Archival ROI -- 4FGL ###############################')
print('######################################################################################')
gta.print_roi()
gta_sources_name=[]
gta_sources_offset=[]
gta_sources_flux_prefit=[]
gta_sources_flux_err_prefit=[]
gta_sources_spec_type_prefit=[]
gta_sources_param_names_prefit=[]
gta_sources_param_val_prefit=[]
gta_sources_param_err_prefit=[]
for i in range(0,len(gta.roi.sources)):
    name=gta.roi.sources[i]['name']
    offset=gta.roi.sources[i]['offset']
    f=gta.roi.sources[i]['flux']
    f_err=gta.roi.sources[i]['flux_err']
    spec_typ=gta.roi.sources[i]['SpectrumType']
    p_name=gta.roi.sources[i]['param_names']
    p_val=gta.roi.sources[i]['param_values']
    p_err=gta.roi.sources[i]['param_errors']
############################################
    gta_sources_name.append(name)
    gta_sources_offset.append(offset)
    gta_sources_flux_prefit.append(f)
    gta_sources_flux_err_prefit.append(f_err)
    gta_sources_spec_type_prefit.append(spec_typ)
    gta_sources_param_names_prefit.append(p_name)
    gta_sources_param_val_prefit.append(p_val)
    gta_sources_param_err_prefit.append(p_err)
############################################
df1=pd.DataFrame({'Name': gta_sources_name,
    'Offset': gta_sources_offset,
    'Flux_pre':gta_sources_flux_prefit,
       'Flux_err_pre':gta_sources_flux_err_prefit,
      'SpectrumType_pre':gta_sources_spec_type_prefit,
       'param_names_pre':gta_sources_param_names_prefit,
        'param_values_pre':gta_sources_param_val_prefit,
       'param_errors_pre':gta_sources_param_err_prefit    
                 })
print('######################################################################################')
print('######################################################################################')
print('#################### Residmap prefit ###############################')
print('######################################################################################')
resid = gta.residmap('residmap_prefit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0}, make_plots=True)

print('######################################################################################')
print('######################################################################################')
print('#################### TSmap prefit ###############################')
print('######################################################################################')
tsmap_prefit = gta.tsmap('tsmap_prefit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0}, make_plots=True)

print('######################################################################################')
print('######################################################################################')
print('#################### Optimize ROI ###############################')
print('######################################################################################')
gta.optimize()

print('######################################################################################')
print('######################################################################################')
print('#################### Optimized ROI ###############################')
print('######################################################################################')
gta.print_roi()


print('######################################################################################')
print('######################################################################################')
print('#################### First source in ROI aka smallest offset ###############################')
print('######################################################################################')
print(gta.roi.sources[0].name)

gta_sources_name=[]
gta_sources_offset=[]
gta_sources_flux_postfit=[]
gta_sources_flux_err_postfit=[]
gta_sources_spec_type_postfit=[]
gta_sources_param_names_postfit=[]
gta_sources_param_val_postfit=[]
gta_sources_param_err_postfit=[]
for i in range(0,len(gta.roi.sources)):
    name=gta.roi.sources[i]['name']
    offset=gta.roi.sources[i]['offset']
    f=gta.roi.sources[i]['flux']
    f_err=gta.roi.sources[i]['flux_err']
    spec_typ=gta.roi.sources[i]['SpectrumType']
    p_name=gta.roi.sources[i]['param_names']
    p_val=gta.roi.sources[i]['param_values']
    p_err=gta.roi.sources[i]['param_errors']
############################################
    gta_sources_name.append(name)
    gta_sources_offset.append(offset)
    gta_sources_flux_postfit.append(f)
    gta_sources_flux_err_postfit.append(f_err)
    gta_sources_spec_type_postfit.append(spec_typ)
    gta_sources_param_names_postfit.append(p_name)
    gta_sources_param_val_postfit.append(p_val)
    gta_sources_param_err_postfit.append(p_err)
############################################
df2=pd.DataFrame({'Name': gta_sources_name,
    'Offset': gta_sources_offset,
    'Flux_post':gta_sources_flux_postfit,
       'Flux_err_post':gta_sources_flux_err_postfit,
      'SpectrumType_post':gta_sources_spec_type_postfit,
       'param_names_post':gta_sources_param_names_postfit,
        'param_values_post':gta_sources_param_val_postfit,
       'param_errors_post':gta_sources_param_err_postfit    
                 })


concentated = pd.merge(df1, df2, on='Name')
concentated['SpectrumType_pre'].values==concentated['SpectrumType_post'].values

delta_flux=[]
delta_index=[]
for i in range(0,len(concentated)):
    r_f=concentated['Flux_post'].values[i]/concentated['Flux_pre'].values[i]
    r_i=concentated['param_values_post'].values[i][1]/concentated['param_values_pre'].values[i][1]
    delta_flux.append(r_f)
    delta_index.append(r_i)

concentated["Delta_Flux"] = [item for item in delta_flux]
concentated["Delta_Index"] = [item for item in delta_index]
outname = 'AAA.csv'
concentated.to_csv(outname, sep="\t", index = False, header=True)
print(concentated)
print('######################################################################################')
print('######################################################################################')
print('#################### Delta (pre- & post-fit) ###############################')
print('######################################################################################')
concentated_deltaflux=concentated.sort_values(by=['Delta_Flux'], ascending=False)
outname = 'Flux&Index_ratios.csv'
concentated_deltaflux.to_csv(outname, sep="\t", index = False, header=True)
print('######################################################################################')
print('######################################################################################')
print('#################### lightcurve ###############################')
print('######################################################################################')
print('Checking whether lightcurve is desired:')
#concentated=concentated.sort_values(by=['ts'], ascending=False)
print(concentated)
print('args.lightcurve = '+str(args.lightcurve))
if args.lightcurve == 'yes':
    which_source=0
    lc = gta.lightcurve(concentated["Name"].values[which_source], nbins=7)
    print((lc['tmin'][0]),(lc['tmax'][0]))
    times = [(lc['tmin_mjd'][0]),(lc['tmax_mjd'][0])]
    t = Time(times, format='mjd', scale='utc')
    t.format = 'isot'
    print(t)
    ###############################################################################################################
    ###############################################################################################################
    fig_sens, ax = plt.subplots(figsize=(10,10))
    #plt.rc('grid', linestyle="dashed", color='grey')
    binsize=1
    ###############################################################################################################
    ###############################################################################################################
    ax1=plt.subplot(211)
    plt.axvline(tmin ,color='green', linestyle='--')
    plt.axvline(tmax ,color='green', linestyle='--')
    tmean = (lc['tmin'] + lc['tmax'])/2
    plt.errorbar(tmean, lc['flux100'],
                 xerr=binsize*24*60*60/2, yerr=lc['flux100_err'], 
                 linestyle='None', marker='o',
                 label='Fermi-LAT > 0.1 GeV '+
                 '\n->'+str(concentated["Name"].values[which_source]))
    plt.ylabel('Flux[cm$^{-2}$ s$^{-1}$]', fontsize=20)
    plt.xlabel('Fermi seconds since 2001.0 UTC')
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.axhline(concentated["Flux_pre"].values[which_source],
                color='k', linestyle='--', label = '4FGL value')
    plt.legend(loc='upper left', fancybox=True, shadow=True)
    plt.grid(True)
    ###############################################################################################################
    ###############################################################################################################
    ax2 = plt.subplot(212, sharex=ax1)
    plt.axvline(tmin ,color='green', linestyle='--')
    plt.axvline(tmax ,color='green', linestyle='--')
    plt.errorbar(tmean, lc['param_values'][:, 1],
                 xerr=binsize*24*60*60/2, yerr=lc['param_errors'][:, 1],
                 linestyle='None', marker='o', label='Fermi-LAT Index '+
                 '\n->'+str(concentated["Name"].values[which_source]))
    plt.ylabel('Photon Index', fontsize=20)
    plt.xlabel('Fermi seconds since 2001.0 UTC')
    plt.axhline(concentated['param_values_pre'].values[which_source][1],
                color='k', linestyle='--', label = '4FGL value')
    plt.legend(loc='lower left', fancybox=True, shadow=True)
    plt.grid(True)
    plt.subplots_adjust(left=0.125,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.01, 
                        hspace=0.01)
    ###############################################################################################################
    ###############################################################################################################
    outname = 'LC_for '+str(concentated["Name"].values[which_source])+'.pdf' 
    plt.savefig(outname)
else:
    print("!!!!!!!!!!->lightcurve not desired<-!!!!!!!!!!")
print('######################################################################################')
print('######################################################################################')
print('#################### Residmap postfit ###############################')
print('######################################################################################')
resid = gta.residmap('residmap_postfit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0}, make_plots=True)
print('######################################################################################')

print('######################################################################################')
print('######################################################################################')
print('#################### TSmap postfit ###############################')
print('######################################################################################')
tsmap_postfit = gta.tsmap('tsmap_postfit',model={'SpatialModel' : 'PointSource', 'Index' : 2.0}, make_plots=True)

print('######################################################################################')
print('######################################################################################')
print('#################### Source Finder ###############################')
print('iterative source-finding algorithm that uses peak detection on a TS map to find new source candidates.')
src = gta.find_sources()
print('######################################################################################')
print('######################################################################################')
print('Done and '+'Total Run-time for fermipy_runner.py is: '+str(time.time() - start_fermipy))
print('######################################################################################')
print('######################################################################################')


