
echo "Running the Binary Experiments (Basic)"
cd ..
nitrs=1
for itr in $(seq 1 $nitrs); do
    echo $itr
    echo "bssize: 16\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_16_itr_$itr"	
    ./docker/run-experiment -p m:17 r17:17 c:17 p:64 >> "repro/data/log_basic_bs_preproc_16_itr_$itr"
    echo "(basic preproc BS size 16), iteration $itr done"
    echo "bssize: 16\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_16_itr_$itr"	
    ./docker/run-experiment  -t 64 bbsearch 16 >> "repro/data/log_basic_bs_online_16_itr_$itr"
    echo "(basic online BS size 16), iteration $itr done"
done    
  
nitrs=1
for itr in $(seq 1 $nitrs); do
    echo "bssize: 18\nis_optimized: 0\n" > 		            "repro/data/log_basic_bs_preproc_18_itr_$itr"	
    ./docker/run-experiment -p m:19 r19:19 c:19 p:64 >> "repro/data/log_basic_bs_preproc_18_itr_$itr"
    echo "(basic preproc BS size 18), iteration $itr done"
    echo "bssize: 18\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_18_itr_$itr"
    ./docker/run-experiment  -t 64 bbsearch 18 >> "repro/data/log_basic_bs_online_18_itr_$itr"
    echo "(basic online BS size 18), iteration $itr done"	
done

nitrs=1
for itr in $(seq 1 $nitrs); do
    echo "bssize: 20\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_20_itr_$itr"	
    ./docker/run-experiment -p m:21 r21:21 c:21 p:64 >> "repro/data/log_basic_bs_preproc_20_itr_$itr"
    echo "(basic preproc BS size 20), iteration $itr done"
    echo "bssize: 20\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_20_itr_$itr"
    ./docker/run-experiment  -t 64 bbsearch 20 >> "repro/data/log_basic_bs_online_20_itr_$itr"
    echo "(basic online BS size 20), iteration $itr done"
done

nitrs=1
for itr in $(seq 1 $nitrs); do
    echo "bssize: 22\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_22_itr_$itr"	
   ./docker/run-experiment -p m:23 r23:23 c:23 p:64 >> "repro/data/log_basic_bs_preproc_22_itr_$itr"
   echo "(basic preproc BS size 22), iteration $itr done"
   echo "bssize: 22\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_22_itr_$itr"
   ./docker/run-experiment  -t 64 bbsearch 22 >>  "repro/data/log_basic_bs_online_22_itr_$itr"
   echo "(basic online BS size 22), iteration $itr done"
done

nitrs=1
for itr in $(seq 1 $nitrs); do
  echo "bssize: 24\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_24_itr_$itr"
  ./docker/run-experiment -p m:25 r25:25 c:25 p:64 >>  "repro/data/log_basic_bs_preproc_24_itr_$itr"
  echo "(basic preproc BS size 24), iteration $itr done"
  echo "bssize: 24\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_24_itr_$itr"
  ./docker/run-experiment  -t 64 bbsearch 24 >>   "repro/data/log_basic_bs_online_24_itr_$itr"
  echo "(basic online BS size 24), iteration $itr done"
done

nitrs=1
for itr in $(seq 1 $nitrs); do
  echo "bssize: 26\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_26_itr_$itr"
  ./docker/run-experiment -p m:27 r27:27 c:27 p:64 >>  "repro/data/log_basic_bs_preproc_26_itr_$itr"
  echo "(basic preproc BS size 26), iteration $itr done"	
  echo "bssize: 26\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_26_itr_$itr"
  ./docker/run-experiment  -t 64 bbsearch 26 >> "repro/data/log_basic_bs_online_26_itr_$itr"
  echo "(online preproc BS size 26), iteration $itr done"
done


echo "Running the Binary Search Experiments (Optimized)"
nitrs=1
for itr in $(seq 1 $nitrs); do
   echo "bssize: 16\nis_optimized: 1\n" 			  >  "repro/data/log_opt_bs_preproc_16_itr_$itr"	 
  ./docker/run-experiment -p  i16:1 c:17 p:64 	  >> "repro/data/log_opt_bs_preproc_16_itr_$itr"
  echo "(optimized preproc BS size 16), iteration $itr done"
   echo "bssize: 16\nis_optimized: 1\n" 				     > "repro/data/log_opt_bs_online_16_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 16   >> "repro/data/log_opt_bs_online_16_itr_$itr"
   echo "(optimized online BS size 16), iteration $itr done"

   echo "bssize: 18\nis_optimized: 1\n" 		     	              > "repro/data/log_opt_bs_preproc_18_itr_$itr"	
  ./docker/run-experiment -p  i18:1 c:19 p:64                        >> "repro/data/log_opt_bs_preproc_18_itr_$itr"
    echo "(optimized preproc BS size 18), iteration $itr done"
   echo "bssize: 18\nis_optimized: 1\n"                                    > "repro/data/log_opt_bs_online_18_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 18   >> "repro/data/log_opt_bs_online_18_itr_$itr"
   echo "(optimized online BS size 18), iteration $itr done"

  echo "bssize: 20\nis_optimized: 1\n"                         > "repro/data/log_opt_bs_preproc_20_itr_$itr"
  ./docker/run-experiment -p  i20:1 c:21 p:64    >> "repro/data/log_opt_bs_preproc_20_itr_$itr"
  echo "(optimized preproc BS size 20), iteration $itr done"
   echo "bssize: 20\nis_optimized: 1\n"                                   > "repro/data/log_opt_bs_online_20_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 20  >> "repro/data/log_opt_bs_online_20_itr_$itr"
  echo "(optimized online BS size 20), iteration $itr done"

   echo "bssize: 22\nis_optimized: 1\n"        > "repro/data/log_opt_bs_preproc_22_itr_$itr"
 ./docker/run-experiment -p  i22:1 c:23 p:64    >> "repro/data/log_opt_bs_preproc_22_itr_$itr"
  echo "(optimized preproc BS size 22), iteration $itr done"
   echo "bssize: 22\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_22_itr_$itr"
 ./docker/run-experiment -t 64 bsearch 22  >> "repro/data/log_opt_bs_online_22_itr_$itr"
  echo "(optimized online BS size 22), iteration $itr done"

 echo "bssize: 24\nis_optimized: 1\n" > "repro/data/log_opt_bs_preproc_24_itr_$itr"
 ./docker/run-experiment -p  i24:1 c:25 p:64 >> "repro/data/log_opt_bs_preproc_24_itr_$itr"
  echo "(optimized preproc BS size 24), iteration $itr done"
   echo "bssize: 24\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_24_itr_$itr"
 ./docker/run-experiment -t 64 bsearch 24  >> "repro/data/log_opt_bs_online_24_itr_$itr"
  echo "(optimized online BS size 24), iteration $itr done"

  echo "bssize: 26\nis_optimized: 1\n" > "repro/data/log_opt_bs_preproc_26_itr_$itr"
 ./docker/run-experiment -p  i26:1 c:27 p:64 >> "repro/data/log_opt_bs_preproc_26_itr_$itr"
  echo "(optimized preproc BS size 26), iteration $itr done"
  echo "bssize: 26\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_26_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 26  >> "repro/data/log_opt_bs_online_26_itr_$itr"
  echo "(optimized online BS size 26), iteration $itr done"

    echo "bssize: 28\nis_optimized: 1\n" > "repro/data/log_opt_bs_preproc_28_itr_$itr"
 ./docker/run-experiment -p  i28:1 c:29 p:64 >> "repro/data/log_opt_bs_preproc_28_itr_$itr"
  echo "(optimized preproc BS size 28), iteration $itr done"
  echo "bssize: 28\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_28_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 28  >> "repro/data/log_opt_bs_online_28_itr_$itr"
  echo "(optimized online BS size 28), iteration $itr done"
done
