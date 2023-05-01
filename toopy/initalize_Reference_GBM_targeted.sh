echo 'Give it a couple of seconds to wake up....;)'
python crossmatch_ranked.py \
-observatory 'Roque de los Muchachos' \
-zenith 50 \
-moon_separation 30 \
-time_res 1 \
-flavour GBM \
-vol 0.9 \
-url 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2022/bn221201517/quicklook/glg_healpix_all_bn221201517.fit' \
-catalog 'IX/67/4fgldr3' \
-ranking_method 'Xmatch' \
-TransNum_TrigID 691590290 \
-too_span 'daily' &