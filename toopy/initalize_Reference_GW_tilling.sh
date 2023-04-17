echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 2 \
-flavour GW \
-vol 0.9 \
-url 'https://gracedb.ligo.org/api/superevents/MS230117q/files/bayestar.multiorder.fits,1' \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'tiled_GW' \
-GraceID 'MS230117q' \
-Rev '1' \
-mode 'diagnostic' \
-instrument_FOV '25' &