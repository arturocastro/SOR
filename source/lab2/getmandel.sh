grep -A 3 "Timer 'Series+monte'" $1 | grep 'Time Per Call' | awk '{ print $NF }'
