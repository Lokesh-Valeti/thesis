
echo "Running the Binary Experiments (Basic)"
cd ..
nitrs=0
for itr in $(seq 1 $nitrs); do
    echo $itr
    echo "bssize: 16\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_16_itr_$itr"	
    ./docker/run-experiment -p m:17 r17:17 c:17 p:64 >> "repro/data/log_basic_bs_preproc_16_itr_$itr"
    echo "preprocessing_bs_16 (basic preproc) - iteration $itr"
    echo "bssize: 16\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_16_itr_$itr"	
    ./docker/run-experiment  -t 64 bbsearch 16 >> "repro/data/log_basic_bs_online_16_itr_$itr"
    echo "preprocessing_bs_16 (basic online) - iteration $itr"
done    
  
nitrs=0
for itr in $(seq 1 $nitrs); do
    echo "bssize: 18\nis_optimized: 0\n" > 		            "repro/data/log_basic_bs_preproc_18_itr_$itr"	
    ./docker/run-experiment -p m:19 r19:19 c:19 p:64 >> "repro/data/log_basic_bs_preproc_18_itr_$itr"
    echo "preprocessing_bs_18 (basic online)"
    echo "bssize: 18\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_18_itr_$itr"
    ./docker/run-experiment  -t 64 bbsearch 16 >> "repro/data/log_basic_bs_online_18_itr_$itr"
    echo "preprocessing_bs_18 (basic online)"	
done

nitrs=0
for itr in $(seq 1 $nitrs); do
    echo "bssize: 20\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_20_itr_$itr"	
    ./docker/run-experiment -p m:21 r21:21 c:21 p:64 >> "repro/data/log_basic_bs_preproc_20_itr_$itr"
    echo "preprocessing_bs_20 (basic preproc)"
    echo "bssize: 20\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_20_itr_$itr"
    ./docker/run-experiment  -t 64 bbsearch 16 >> "repro/data/log_basic_bs_online_20_itr_$itr"
    echo "preprocessing_bs_20 (basic online)"
done

nitrs=0
for itr in $(seq 1 $nitrs); do
    echo "bssize: 22\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_22_itr_$itr"	
   ./docker/run-experiment -p m:23 r23:23 c:23 p:64 >> "repro/data/log_basic_bs_preproc_22_itr_$itr"
   echo "preprocessing_bs_22 (basic preproc)"
   echo "bssize: 22\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_22_itr_$itr"
   ./docker/run-experiment  -t 64 bbsearch 16 >>  "repro/data/log_basic_bs_online_22_itr_$itr"
   echo "preprocessing_bs_22 (basic online)"
done

nitrs=0
for itr in $(seq 1 $nitrs); do
  echo "bssize: 24\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_24_itr_$itr"
  ./docker/run-experiment -p m:25 r25:25 c:25 p:64 >>  "repro/data/log_basic_bs_preproc_24_itr_$itr"
  echo "preprocessing_bs_24 (basic preproc)"
  echo "bssize: 24\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_24_itr_$itr"
  ./docker/run-experiment  -t 64 bbsearch 16 >>   "repro/data/log_basic_bs_online_24_itr_$itr"
  echo "preprocessing_bs_24 (basic online)"
done

nitrs=0
for itr in $(seq 1 $nitrs); do
  echo "bssize: 26\nis_optimized: 0\n" > "repro/data/log_basic_bs_preproc_26_itr_$itr"
  ./docker/run-experiment -p m:27 r27:27 c:27 p:64 >>  "repro/data/log_basic_bs_preproc_26_itr_$itr"
  echo "bssize: 26\nis_optimized: 0\n" > "repro/data/log_basic_bs_online_26_itr_$itr"
  ./docker/run-experiment  -t 64 bbsearch 16 >> "repro/data/log_basic_bs_online_26_itr_$itr"
  echo "preprocessing_bs_26 (basic online)"
done


echo "Running the Heap Extract Experiments (Optimized)"
nitrs=2
for itr in $(seq 1 $nitrs); do
   echo "bssize: 16\nis_optimized: 1\n" 			  >  "repro/data/log_opt_bs_preproc_16_itr_$itr"	 
  ./docker/run-experiment -p  i16:1 c:17 p:64 	  >> "repro/data/log_opt_bs_preproc_16_itr_$itr"
  echo "preprocessing_bs_16 (opt preproc)"
   echo "bssize: 16\nis_optimized: 1\n" 				     > "repro/data/log_opt_bs_online_16_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 16   >> "repro/data/log_opt_bs_online_16_itr_$itr"
  echo "preprocessing_bs_16 (opt online)"

   echo "bssize: 18\nis_optimized: 1\n" 		     	              > "repro/data/log_opt_bs_preproc_18_itr_$itr"	
  ./docker/run-experiment -p  i18:1 c:19 p:64                        >> "repro/data/log_opt_bs_preproc_18_itr_$itr"
   echo "preprocessing_bs_18 (opt preproc)"
   echo "bssize: 18\nis_optimized: 1\n"                                    > "repro/data/log_opt_bs_online_18_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 18   >> "repro/data/log_opt_bs_online_18_itr_$itr"
  echo "preprocessing_bs_18 (opt online)"

  echo "bssize: 20\nis_optimized: 1\n"                         > "repro/data/log_opt_bs_preproc_20_itr_$itr"
  ./docker/run-experiment -p  i20:1 c:21 p:64    >> "repro/data/log_opt_bs_preproc_20_itr_$itr"
  echo "preprocessing_bs_20 (opt preproc)"
   echo "bssize: 20\nis_optimized: 1\n"                                   > "repro/data/log_opt_bs_online_20_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 20  >> "repro/data/log_opt_bs_online_20_itr_$itr"
  echo "preprocessing_bs_20 (opt online)"

   echo "bssize: 22\nis_optimized: 1\n"        > "repro/data/log_opt_bs_preproc_22_itr_$itr"
 ./docker/run-experiment -p  i22:1 c:23 p:64    >> "repro/data/log_opt_bs_preproc_22_itr_$itr"
  echo "preprocessing_bs_22 (opt preproc)"
   echo "bssize: 22\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_22_itr_$itr"
 ./docker/run-experiment -t 64 bsearch 22  >> "repro/data/log_opt_bs_online_22_itr_$itr"
  echo "preprocessing_bs_22 (opt online)"

 echo "bssize: 24\nis_optimized: 1\n" > "repro/data/log_opt_bs_preproc_24_itr_$itr"
 ./docker/run-experiment -p  i24:1 c:25 p:64 >> "repro/data/log_opt_bs_preproc_24_itr_$itr"
  echo "preprocessing_bs_24 (opt preproc)"
   echo "bssize: 24\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_24_itr_$itr"
 ./docker/run-experiment -t 64 bsearch 24  >> "repro/data/log_opt_bs_online_24_itr_$itr"
  echo "preprocessing_bs_24 (opt online)"

  echo "bssize: 26\nis_optimized: 1\n" > "repro/data/log_opt_bs_preproc_26_itr_$itr"
 ./docker/run-experiment -p  i26:1 c:27 p:64 >> "repro/data/log_opt_bs_preproc_26_itr_$itr"
  echo "preprocessing_bs_26 (opt preproc)"
  echo "bssize: 26\nis_optimized: 1\n" > "repro/data/log_opt_bs_online_26_itr_$itr"
  ./docker/run-experiment -t 64 bsearch 26  >> "repro/data/log_opt_bs_online_26_itr_$itr"
  echo "preprocessing_bs_26 (opt online)"
done
