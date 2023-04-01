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
-ranking_method 'fermipy' \
-fermitools_refdata_path '/opt/anaconda3/envs/toopy/share/fermitools/refdata' \
-Rev IC170922A \
-lightcurve 'no' &
