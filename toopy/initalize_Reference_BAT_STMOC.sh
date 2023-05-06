echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour BAT \
-ra 266.933 \
-dec '-68.260' \
-error 0.1 \
-obs_night 59914.51831019 \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'STMOC' \
-TransNum_TrigID 1142847 \
-too_span 'daily' &
