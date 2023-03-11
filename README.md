# ToOpy
# -------------------------
# -------------------------
## Setting up the environment
### fermipy compatible
Setting up conda environment with the right dependencies:

$ `git clone https://github.com/arterom/ToOpy.git`

$ `git clone https://github.com/fermiPy/fermipy.git`

$ `cd fermipy`

$ `cp ../toopy/conda_environment.yml .` Or copy "conda_environment.yml" file from toopy folder to newly created fermipy folder.

$ `conda env create -f conda_environment.yml` ***

$ `conda activate toopy`

$ `python setup.py install` Installing fermipy package

# -------------------------
# -------------------------

## Quickstart Guide

$ `cd ../toopy`

$ `conda activate toopy`

## User preferences
Edit lines 23-27 in "ToOpy_gcn_listener_main.py" according to user preferences:
```
######################################################################################
observatory='"Roque de los Muchachos"' #***
max_zenith='50'
moon_separation='30'
time_resolution='1' # time in hours per observing slot
fermitools_refdata_path='/opt/anaconda3/envs/toopy/share/fermitools/refdata'
######################################################################################
#***To get all possible sites use astropy-module:
#from astropy.coordinates import EarthLocation
#EarthLocation.get_site_names()
```


## Listeing to alerts
$ `python ToOpy_gcn_listener_main.py`

Notice the constant stream of hourly GW alerts @ x:15 (Revision 0) & x:20 (Revision 1) ;D
Plus, every now and then some FermiGBM alerts ('FERMI_GBM_FIN_POS' & 'FERMI_GBM_SUBTHRESH')


# -------------------------
# -------------------------

## Testing Guide


$ `cd ../toopy`

$ `conda activate toopy`

### Option 1; Testing pipeline with test alerts from ".xml" files
Step 1: Modifiy "ToOpy_gcn_listener_main.py"

-> Comment line 324

-> Uncomment a given alert (lines 324-333)

-> Uncomment line 334 & 335

```
# Listen for VOEvents until killed with Control-C.
#gcn.listen(handler=handler)

# Templates for test alerts
#payload = open('./xml_test_alerts/FermiGBM/test_GBM_SubTresh.xml', 'rb').read()
#payload = open('./xml_test_alerts/FermiGBM/test_GBM_Final.xml', 'rb').read()
#payload = open('./xml_test_alerts/GW/test_GW_nuniq.xml', 'rb').read()
#payload = open('./xml_test_alerts/GW/test_GW_nside.xml', 'rb').read()
#payload = open('./xml_test_alerts/IceCube_TRACK/GOLD_Aug11_Rev1.xml', 'rb').read()
payload = open('./xml_test_alerts/IceCube_CASCADE/CASCADE_Apr.xml', 'rb').read()

root = lxml.etree.fromstring(payload)
handler(payload, root)
```

Step 2: $ `python ToOpy_gcn_listener_main.py`

### Option 2; Testing pipeline with test alerts from script (for Track-like alerts)
Step 1: Edit "initalize_custom_IceCube_TRACK_alert.sh" according to user preferences:
```
echo 'Give it a second or two to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour Track \
-ra 266.80 \
-dec '-3.58' \
-error 1 \
-obs_night 59620 \
-catalog 'J/ApJS/247/33/4fgl' \
-ranking_method 'VarInd' \
-fermitools_refdata_path '/opt/anaconda3/envs/mma_broker/share/fermitools/refdata' \
-Rev manually &
```

Step 2: $ `./initalize_custom_IceCube_TRACK_alert.sh`

# -------------------------
# -------------------------




