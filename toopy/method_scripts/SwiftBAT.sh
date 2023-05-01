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
printf 'Running Xmatch with method: %s\n' "$9"
python crossmatch_ranked.py \
-observatory "$1" \
-zenith "$2" \
-moon_separation "$3" \
-time_res "$4" \
-flavour BAT \
-ra "$5" \
-dec " $6" \
-error "$7" \
-obs_night "$8" \
-catalog 'IX/67/4fgldr3' \
-ranking_method "$9" \
-TransNum_TrigID "${10}" \
-too_span "${11}" &
