input_bench="ispd18_sample"
bench_dir="../../../ispd18"
tmp="ispd18_test1_net1_net254"
n=1
while [ $n -le 1 ]
do
    echo "Welcome $n times."
    for i in $input_bench
    do 
        echo "$i ..."
        rm $bench_dir/$i/$i\_log\_mov\_dr.txt
        echo "detailed router"
        ./TritonRoute -lef $bench_dir/$i/$i.input.lef -def $bench_dir/$i/$i.input.def -guide $bench_dir/$i/$i.input.guide -output $bench_dir/$i/$i.out.dr.def -threads 8  >> $bench_dir/$i/$i\_log\_dr.txt
        echo "$i done!"
    done 
    n=$(( n+1 ))	 # increments $n
done
