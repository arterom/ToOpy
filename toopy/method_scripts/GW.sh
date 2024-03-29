printf 'Arg 1 is: %s\n' "$1"
printf 'Arg 2 is: %s\n' "$2"
printf 'Arg 3 is: %s\n' "$3"
printf 'Arg 4 is: %s\n' "$4"
printf 'Arg 5 is: %s\n' "$5"
printf 'Arg 6 is: %s\n' "$6"
printf 'Arg 7 is: %s\n' "$7"
printf 'Arg 8 is: %s\n' "$8"
printf 'Arg 9 is: %s\n' "$9"
printf 'Arg 10 is: %s\n' "${10}"
printf 'Arg 11 is: %s\n' "${11}"
printf 'Running Xmatch with method: %s\n' "$6"
python crossmatch_ranked.py \
-observatory "$1" \
-zenith "$2" \
-moon_separation "$3" \
-time_res "$4" \
-flavour GW \
-vol 0.9 \
-url "$5" \
-catalog 'IX/67/4fgldr3' \
-ranking_method "$6" \
-GraceID "$7" \
-Rev "$8" \
-mode "$9" \
-instrument_FOV "${10}" \
-too_span "${11}" &
