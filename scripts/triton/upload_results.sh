# costs:        wl_cost = 0;
#         cong_cost = 0;
#         via_cost = 0;
#         via_abs_cost = 0;
#         wl_abs_cost = 0;
#         path_cost = 0;
#         alignment_cost = 0;
#         hpwl_cost = 0;
# input_bench="ispd18_test1 ispd18_test2 \
#  ispd18_test3 ispd18_test4 \
#  ispd18_test5 ispd18_test6 \
#  ispd18_test7 ispd18_test8 \
#  ispd18_test9 ispd18_test10 \
#  ispd19_test1 ispd19_test2 \
#  ispd19_test3 ispd19_test4 \
#  ispd19_test5 ispd19_test6 \
#  ispd19_test7 ispd19_test8 \
#  ispd19_test9 ispd19_test10"
# -def $bench_dir/$i/$i.input.cadence.def \

#  input_bench="ispd18_sample ispd18_sample2 ispd18_sample3 \
#   ispd18_test1 ispd18_test2 \
#   ispd18_test3 ispd18_test4 \
#   ispd18_test5 ispd18_test6 \
#   ispd18_test7 ispd18_test8 \
#   ispd18_test9 ispd18_test10"

#  input_bench="ispd18_test1 ispd18_test2 \
#  ispd18_test3 ispd18_test4 \
#  ispd18_test5 ispd18_test6 \
#  ispd18_test7 ispd18_test8 \
#  ispd18_test9 ispd18_test10"

#input_bench="ispd18_sample ispd18_sample2 ispd18_sample3"
# -def $bench_dir/$i/$i.input.cadence.inst4183.def \
input_bench="ispd18_sample"
bench_dir="../../../ispd18"
tmp="ispd18_test1_net1_net254"
experiment="15_experiment"
mkdir "../reports/$experiment"
cp run.sh ../reports/$experiment
n=1
while [ $n -le 1 ]
do
    echo "Welcome $n times."
    for i in $input_bench
    do 
        echo "$i ..."
        cp $bench_dir/$i/$i.out.dr.def $bench_dir/$i/$i.out.mov.dr.def
        cp $bench_dir/$i/$i.input.guide $bench_dir/$i/$i.solution.guide	 
        echo "upload results to server"
        ./copy_to_cadence_server.sh $bench_dir/$i/$i.out.mov.dr.def ispd18eval/benchmarks/$i/
        ./copy_to_cadence_server.sh $bench_dir/$i/$i.solution.guide ispd18eval/benchmarks/$i/

         echo "end copy results"

         echo "$i done!"
    done 
    n=$(( n+1 ))	 # increments $n
done
# for i in $input_bench
# do 
#     echo "$i ..."
#     rm $bench_dir/$tmp/$i\_log\_mov.txt
#     # ./iccad19gr -lef $bench_dir/$i/$i.input.lef -def $bench_dir/$i/$i.input.def -output $bench_dir/$i/$i.solution.guide -threads 8 >> $bench_dir/$i/$i\_log\_mov.txt
#     ./iccad19gr -lef $bench_dir/$tmp/$i.input.lef -def $bench_dir/$tmp/$i.input.def -output $bench_dir/$tmp/$i.solution.guide -threads 1 >> $bench_dir/$tmp/$i\_log\_mov.txt
#     # ./iccad19gr -lef $bench_dir/$tmp/$i.input.lef -def $bench_dir/$tmp/out.def -output $bench_dir/$tmp/$i.solution.guide -threads 8 >> $bench_dir/$tmp/$i\_log\_mov.txt
#     echo "$i done!"

# done 

