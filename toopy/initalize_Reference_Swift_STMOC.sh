echo 'Give it a couple of seconds to wake up....;)'
#from astropy.time import Time
#from astropy.time import TimeDelta

#Event_MJD=59914
#Event_SOD=44705.68

#start_time_mjd = Time(Event_MJD, format='mjd', scale='utc')
#dt_sod = TimeDelta(Event_SOD, format='sec')
#start_time_mjd_actual = start_time_mjd + dt_sod'''

python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour BAT \
-ra 266.933 \
-dec '-68.260' \
-error 3 \
-obs_night 59914.51831019 \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'Xmatch' \
-TransNum_TrigID 1142847 &