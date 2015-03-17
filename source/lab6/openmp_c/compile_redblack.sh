#echo $1;

sizes=( 20 50 100 200 300 400 500 600 750)

for s in "${sizes[@]}"
do
    #echo $s;
    make clean -s;
    make SIZE=$s sorrb;
done
