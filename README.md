# ToOpy
The algorithm is aimed at combining the functionality of a data Broker and observation Scheduler. The Broker intercepts alert notices that are injected into the Gamma Coordination Network (GCN) stream by a range of observatories, and enrich therein content with additional information from galaxy catalogs as well as contemporaneous satelite-data products, e.g. Swift, Fermi-LAT & Fermi-GBM coverage. Once promising astrophysical objects and/or sky regions are ranked in order of priority, the Broker provides the Scheduler with a list of galaxies or sky tiles, that based on user preferences are converted into pointing sequences for automatised Target of Opportunity (ToO) observations.
# -------------------------
# -------------------------
## Setting up the environment
Setting up conda environment with the right dependencies:

$ `git clone https://github.com/arterom/ToOpy.git`

$ `git clone https://github.com/fermiPy/fermipy.git`

$ `cd fermipy`

$ `cp ../toopy/conda_environment.yml .`

$ `conda env create -f conda_environment.yml`

$ `conda activate toopy`

$ `python setup.py install` ## installing fermipy package

# -------------------------
# -------------------------

## Quickstart Guide

$ `cd ../toopy`

$ `conda activate toopy`

## User preferences
Edit lines 27-43 in "toopy_gcn_listener_main.py" according to user preferences:
```
# General Constraints
#############################################
#############################################
observatory='"Roque de los Muchachos"' #***
max_zenith='50'
moon_separation='30'
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

```


## Listeing to alerts
$ `python toopy_gcn_listener_main.py`

Notice the constant stream of hourly GW alerts @ x:15 (Revision 0) & x:20 (Revision 1) ;D
Plus, every now and then some FermiGBM alerts ('FERMI_GBM_FIN_POS' & 'FERMI_GBM_SUBTHRESH')


# -------------------------
# -------------------------

## Testing Guide


$ `cd ../toopy`

$ `conda activate toopy`

### Option 1; Testing pipeline with test alerts from ".xml" files
Step 1: Modifiy "toopy_gcn_listener_main.py"

-> Comment line '444'

-> Uncomment a given alert (lines 447-450)

-> Uncomment line 453 & 454

```
# Listen for VOEvents until killed with Control-C.
#gcn.listen(handler=handler)

# Templates for test alerts
#payload = open('./xml_test_alerts/IceCube_CASCADE/2021/CASCADE_Dez21.xml', 'rb').read()

root = lxml.etree.fromstring(payload)
handler(payload, root)
```

Step 2: $ `python toopy_gcn_listener_main.py`

### Option 2; Testing pipeline with test alerts from script (for Track-like alerts)
Step 1: Edit "initalize_Reference_IC170922A.sh" according to user preferences:
```
echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour Track \
-ra 77.2853 \
-dec '+5.7517' \
-error 0.416388883335 \
-obs_night 58018 \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'VarInd' \
-fermitools_refdata_path '/opt/anaconda3/envs/toopy/share/fermitools/refdata' \
-Rev IC170922A \
-lightcurve 'no' &
```

Step 2: $ `./initalize_Reference_IC170922A.sh`

# -------------------------
# -------------------------




