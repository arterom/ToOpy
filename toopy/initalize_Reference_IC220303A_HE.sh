echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour Track \
-ra 267.8000 \
-dec '+11.4199' \
-error 1.1798333333333335 \
-obs_night 59641 \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'fermipy' \
-fermitools_refdata_path '/opt/anaconda3/envs/toopy/share/fermitools/refdata' \
-Rev IC220303A_1  \
-lightcurve 'yes' \
-too_span 'daily' &
