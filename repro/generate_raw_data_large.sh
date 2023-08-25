
echo "Running the Heap Extract Experiments (Basic)"
cd ..
nitrs=2
for itr in $(seq 1 $nitrs); do
    echo $itr
    echo "heapsize: 28\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_preproc_28_itr_$itr"	
    ./docker/run-experiment -p m:111 a:28 s:27 r28:12 c:56 p:128 >> "repro/data/log_basic_heap_extract_preproc_28_itr_$itr"
    
    for npre in $(seq 1 15); do
	 echo "heapsize: 28\nis_optimized: 0\n" > "repro/data/log_"$npre"_basic_heap_extract_preproc_28_itr_$itr"
    	 ./docker/run-experiment -a -p r28:10 >> "repro/data/log_"$npre"_basic_heap_extract_preproc_28_itr_$itr"
     done 
    echo "preprocessing_heap_28 (basic preproc) - iteration $itr"
    echo "heapsize: 28\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_28_itr_$itr"	
    ./docker/run-experiment heap -m 28 -d 28 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_28_itr_$itr"
    echo "preprocessing_heap_16 (basic online) - iteration $itr"
done    
 
nitrs=2
for itr in $(seq 1 $nitrs); do
    echo "heapsize: 30\nis_optimized: 0\n" > 		            "repro/data/log_basic_heap_extract_preproc_30_itr_$itr"	
    ./docker/run-experiment -p m:119 a:30 s:29 r30:17 c:60 p:128 >> "repro/data/log_basic_heap_extract_preproc_30_itr_$itr"
    echo "preprocessing_heap_30 (basic online)"
    for npre in $(seq 1 15); do
    echo "heapsize: 30\nis_optimized: 0\n" >                        "repro/data/log3_basic_heap_extract_preproc_30_itr_$itr"
    ./docker/run-experiment -a -p r30:10  >> "repro/data/log_"$npre"_basic_heap_extract_preproc_30_itr_$itr"
     echo "preprocessing_heap_30 (basic online)"
	done
    echo "heapsize: 30\nis_optimized: 0\n" > "repro/data/log_basic_heap_extract_online_30_itr_$itr"
    ./docker/run-experiment heap -m 30 -d 30 -i 0 -e 1 -opt 0 -s 0 -itr itr >> "repro/data/log_basic_heap_extract_online_30_itr_$itr"
    echo "preprocessing_heap_30 (basic online)"	
done


echo "Running the Heap Extract Experiments (Optimized)"
nitrs=2
for itr in $(seq 1 $nitrs); do
     echo "heapsize: 28\nis_optimized: 1\n" 				  >  "repro/data/log_opt_heap_extract_preproc_28_itr_$itr"	 
     ./docker/run-experiment -p  m:111 a:28 s:27 i27.3:1 c:56 p:128 	  >> "repro/data/log_opt_heap_extract_preproc_28_itr_$itr"
     echo "preprocessing_heap_16 (opt preproc)"
     echo "heapsize: 28\nis_optimized: 1\n" 				     > "repro/data/log_opt_heap_extract_online_28_itr_$itr"
     ./docker/run-experiment heap  -m 28 -d 28 -i 0 -e 1 -opt 1 -s 0 -itr itr  >> "repro/data/log_opt_heap_extract_online_28_itr_$itr"
     echo "preprocessing_heap_16 (opt online)"

     echo "heapsize: 30\nis_optimized: 1\n" 		     	              > "repro/data/log_opt_heap_extract_preproc_30_itr_$itr"	
     ./docker/run-experiment -p m:119 a:30 s:29 i29.3:1 c:60 p:128               >> "repro/data/log_opt_heap_extract_preproc_30_itr_$itr"
     echo "preprocessing_heap_18 (opt preproc)"
     echo "heapsize: 30\nis_optimized: 1\n"                                    > "repro/data/log_opt_heap_extract_online_30_itr_$itr"
     ./docker/run-experiment heap  -m 30 -d 30 -i 0 -e 1 -opt 1 -s 0 -itr itr  >> "repro/data/log_opt_heap_extract_online_30_itr_$itr"
     echo "preprocessing_heap_30 (opt online)"
  done
