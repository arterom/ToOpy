printf 'Arg 1 is: %s\n' "$1"
printf 'THIS WILL LAUNCH STMOC INTERCEPTOR'
python --version
python stmoc_interceptor.py \
-string "$1" &
