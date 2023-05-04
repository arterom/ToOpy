echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour Track \
-ra 197.448776  \
-dec '-23.383831' \
-error 3 \
-obs_night 57986 \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'fermipy' \
-fermitools_refdata_path '/opt/anaconda3/envs/mma_broker/share/fermitools/refdata' \
-Rev GW170817_GRB170817A \
-lightcurve 'no' &
