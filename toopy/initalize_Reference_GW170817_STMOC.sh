echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 2 \
-flavour GW \
-vol 0.9 \
-url 'https://dcc.ligo.org/public/0146/G1701985/001/bayestar_no_virgo.fits.gz' \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'STMOC' \
-GraceID 'GW170817' \
-Rev '1' \
-mode 'diagnostic' \
-instrument_FOV '25' \
-too_span 'daily' &
