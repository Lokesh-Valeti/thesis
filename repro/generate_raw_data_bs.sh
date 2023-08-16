 


./docker/run-experiment -p m:17 r17:17 c:17 p:64 > preprocessing_binary_basic_16
./docker/run-experiment -t 64 bbsearch 16 > online_binary_basic_16
echo "online_binary_basic_16"

./docker/run-experiment -p m:19 r19:19 c:19 p:64 > preprocessing_binary_basic_18
./docker/run-experiment -t 64 bbsearch 18 > online_binary_basic_18
echo "online_binary_basic_18"

./docker/run-experiment -p m:21 r21:21 c:21 p:64 > no_latency_preprocessing_binary_basic_20
./docker/run-experiment -t 64 bbsearch 20 > no_latency_online_binary_basic_20
echo "online_binary_basic_20"

./docker/run-experiment -p m:25 r23:23 c:23 p:64 > preprocessing_binary_basic_22
./docker/run-experiment -t 64 bbsearch 22 > online_binary_basic_22
 echo "online_binary_basic_22"

./docker/run-experiment -p m:25 r25:25 c:25 p:64 > preprocessing_binary_basic_24
./docker/run-experiment -t 64 bbsearch 24 > online_binary_basic_24
echo "online_binary_basic_24"

./docker/run-experiment -p m:27 r27:27 c:27 p:64 > no_latency_preprocessing_binary_basic_26
./docker/run-experiment -t 64 bbsearch 26 > no_latency_online_binary_basic_26
 echo "online_binary_basic_26"


echo "Running the Binary Experiments (Basic)"
cd ..
nitrs=10
for itr in $(seq 1 $nitrs); do
    echo $itr
    echo "heapsize: 16\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_16_itr_$itr"	
    ./docker/run-experiment -p m:63 a:16 s:15 r16:90 c:32 p:128 >> "repro/data/log_basic_heap_extract_preproc_16_itr_$itr"
    echo "preprocessing_heap_16 (basic preproc) - iteration $itr"
    echo "heapsize: 16\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_16_itr_$itr"	
    ./docker/run-experiment heap -m 16 -d 16 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_16_itr_$itr"
    echo "preprocessing_heap_16 (basic online) - iteration $itr"
done    
  
nitrs=10
for itr in $(seq 1 $nitrs); do
    echo "heapsize: 18\nis_optimized: 0\n" > 		            "repro/data/log_basic_heap_extract_preproc_18_itr_$itr"	
    ./docker/run-experiment -p m:71 a:18 s:17 r18:102 c:36 p:128 >> "repro/data/log_basic_heap_extract_preproc_18_itr_$itr"
    echo "preprocessing_heap_18 (basic online)"
    echo "heapsize: 18\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_18_itr_$itr"
    ./docker/run-experiment heap -m 18 -d 18 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_18_itr_$itr"
    echo "preprocessing_heap_18 (basic online)"	
done

nitrs=10
for itr in $(seq 1 $nitrs); do
    echo "heapsize: 20\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_20_itr_$itr"	
    ./docker/run-experiment -p m:80 a:20 s:19 r20:114 c:40 p:128 >> "repro/data/log_basic_heap_extract_preproc_20_itr_$itr"
    echo "preprocessing_heap_20 (basic preproc)"
    echo "heapsize: 20\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_20_itr_$itr"
    ./docker/run-experiment -t 8 heap -m 20 -d 20 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_20_itr_$itr"
    echo "preprocessing_heap_20 (basic online)"
done

nitrs=10
for itr in $(seq 1 $nitrs); do
    echo "heapsize: 22\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_22_itr_$itr"	
   ./docker/run-experiment -p m:87 a:22 s:21 r22:126 c:44 p:128 >> "repro/data/log_basic_heap_extract_preproc_22_itr_$itr"
   echo "preprocessing_heap_22 (basic preproc)"
   echo "heapsize: 22\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_22_itr_$itr"
   ./docker/run-experiment -t 8 heap -m 22 -d 22 -i 0 -e 1 -opt 0 -s 0 -itr itr >>  "repro/data/log_basic_heap_extract_online_22_itr_$itr"
   echo "preprocessing_heap_22 (basic online)"
done

nitrs=10
for itr in $(seq 1 $nitrs); do
  echo "heapsize: 24\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_24_itr_$itr"
  ./docker/run-experiment -p m:69 a:23 s:22 r24:132 c:46 >>  "repro/data/log_basic_heap_extract_preproc_24_itr_$itr"
  echo "preprocessing_heap_24 (basic preproc)"
  echo "heapsize: 24\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_24_itr_$itr"
  ./docker/run-experiment -t 8 heap -m 24 -d 24 -i 0 -e 1 -opt 0 -s 0 -itr itr >>   "repro/data/log_basic_heap_extract_online_24_itr_$itr"
  echo "preprocessing_heap_24 (basic online)"
done

nitrs=0
for itr in $(seq 1 $nitrs); do
  echo "heapsize: 26\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_26_itr_$itr"
  ./docker/run-experiment -p m:75 a:25 s:24 r26:100 c:50 p:8 >>  "repro/data/log_basic_heap_extract_preproc_26_itr_$itr"
  echo "heapsize: 26\nis_optimized: 0\n" > "repro/data/logb_basic_heap_extract_preproc_26_itr_$itr"
  ./docker/run-experiment -p -a r26:44  p:8 >>  "repro/data/logb_basic_heap_extract_prproc_26_itr_$itr"
  echo "preprocessing_heap_26 (basic preproc)"
  echo "heapsize: 26\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_26_itr_$itr"
  ./docker/run-experiment -t 8 heap -m 26 -d 26 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_26_itr_$itr"
  echo "preprocessing_heap_26 (basic online)"
done


echo "Running the Heap Extract Experiments (Optimized)"

for itr in $(seq 1 $nitrs); do
   echo "heapsize: 16\nis_optimized: 1\n" 				  >  "repro/data/log_opt_heap_extract_preproc_16_itr_$itr"	 
  ./docker/run-experiment -p  m:63 a:16 s:15  i15.3:1 c:32 p:128 	  >> "repro/data/log_opt_heap_extract_preproc_16_itr_$itr"
  echo "preprocessing_heap_16 (opt preproc)"
   echo "heapsize: 16\nis_optimized: 1\n" 				     > "repro/data/log_opt_heap_extract_online_16_itr_$itr"
  ./docker/run-experiment heap  -m 16 -d 16 -i 0 -e 1 -opt 1 -s 0 -itr itr  >> "repro/data/log_opt_heap_extract_online_16_itr_$itr"
  echo "preprocessing_heap_16 (opt online)"

   echo "heapsize: 18\nis_optimized: 1\n" 		     	              > "repro/data/log_opt_heap_extract_preproc_18_itr_$itr"	
  ./docker/run-experiment -p m:71 a:18 s:17 i17.3:1 c:36 p:128               >> "repro/data/log_opt_heap_extract_preproc_18_itr_$itr"
   echo "preprocessing_heap_18 (opt preproc)"
   echo "heapsize: 18\nis_optimized: 1\n"                                    > "repro/data/log_opt_heap_extract_online_18_itr_$itr"
  ./docker/run-experiment heap  -m 18 -d 18 -i 0 -e 1 -opt 1 -s 0 -itr itr  >> "repro/data/log_opt_heap_extract_online_18_itr_$itr"
  echo "preprocessing_heap_18 (opt online)"

 echo "heapsize: 20\nis_optimized: 1\n"                         > "repro/data/log_opt_heap_extract_preproc_20_itr_$itr"
  ./docker/run-experiment -p m:80 a:20 s:19 i19.3:1 c:40 p:128 >> "repro/data/log_opt_heap_extract_preproc_20_itr_$itr"
  echo "preprocessing_heap_20 (opt preproc)"
   echo "heapsize: 20\nis_optimized: 1\n"                                   > "repro/data/log_opt_heap_extract_online_20_itr_$itr"
  ./docker/run-experiment heap  -m 20 -d 20 -i 0 -e 1 -opt 1 -s 0 -itr itr >> "repro/data/log_opt_heap_extract_online_20_itr_$itr"
  echo "preprocessing_heap_20 (opt online)"

   echo "heapsize: 22\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_preproc_22_itr_$itr"
 ./docker/run-experiment -p m:87 a:22 s:21 i21.3:1 c:44 p:128  >> "repro/data/log_opt_heap_extract_preproc_22_itr_$itr"
  echo "preprocessing_heap_22 (opt preproc)"
   echo "heapsize: 22\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_online_22_itr_$itr"
 ./docker/run-experiment -t 16 heap  -m 22 -d 22 -i 0 -e 1 -opt 1 -s 1 -itr itr >> "repro/data/log_opt_heap_extract_online_22_itr_$itr"
  echo "preprocessing_heap_22 (opt online)"

 echo "heapsize: 24\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_preproc_24_itr_$itr"
 ./docker/run-experiment -p m:95 a:24 s:23 i23.3:1 c:48 p:128 >> "repro/data/log_opt_heap_extract_preproc_24_itr_$itr"
  echo "preprocessing_heap_24 (opt preproc)"
   echo "heapsize: 24\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_online_24_itr_$itr"
 ./docker/run-experiment -t 32 heap  -m 24 -d 24 -i 0 -e 1 -opt 1 -s 0 -itr itr >> "repro/data/log_opt_heap_extract_online_24_itr_$itr"
  echo "preprocessing_heap_24 (opt online)"

  echo "heapsize: 26\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_preproc_26_itr_$itr"
 ./docker/run-experiment -p m:103 a:26 s:25 i25.3:1 c:52 p:128 >> "repro/data/log_opt_heap_extract_preproc_26_itr_$itr"
  echo "preprocessing_heap_26 (opt preproc)"
  echo "heapsize: 26\nis_optimized: 1\n" > "repro/data/log_opt_heap_extract_online_26_itr_$itr"
  ./docker/run-experiment -t 64 heap  -m 26 -d 26 -i 0 -e 1 -opt 1 -s 0 -itr itr >> "repro/data/log_opt_heap_extract_online_26_itr_$itr"
  echo "preprocessing_heap_26 (opt online)"
done


