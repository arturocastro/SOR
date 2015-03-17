grep "Total" $1 | awk '{ print $NF }' | sort | head -n1
