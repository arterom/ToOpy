printf 'Arg 1 is: %s\n' "$1"
printf 'Arg 2 is: %s\n' "$2"
python --version
python fermipy_runner.py \
-outpath "$1" \
-lightcurve "$2" &
