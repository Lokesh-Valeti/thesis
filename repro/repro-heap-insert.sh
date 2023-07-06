 ./run-experiment -p m:16 c:16 > basic_heap_insert_preproc_16
 ./run-experiment -t 64 heap -m 17 -d 16 -i 1 -e 0 -opt 0 -s 0 > log_heap_insert_basic_online_16
 echo "basic heap insert for 2^16"
 
 ./run-experiment -p m:18 c:18 > basic_heap_insert_preproc_18
 ./run-experiment -t 64 heap -m 19 -d 18 -i 1 -e 0 -opt 0 -s 0 > log_heap_insert_basic_online_18
 echo "basic heap insert for 2^18"
 
 ./run-experiment -p  m:20 c:20 > log_heap_insert_basic_pre_20
 ./run-experiment heap -m 21 -d 20 -i 1 -e 0 -opt 0 -s 0 > log_heap_insert_basic_online_20
 echo "basic heap insert for 2^20"

 
 ./run-experiment -p m:22 c:22 > log_heap_insert_basic_pre_22
 ./run-experiment  heap -m 23 -d 22 -i 1 -e 0 -opt 0 -s 0 > log_heap_insert_basic_online_22
 echo "basic heap insert for 2^22"

 ./run-experiment -p m:24 c:24 > log_heap_insert_basic_pre_24
 ./run-experiment heap -m 25 -d 24 -i 1 -e 0 -opt 0 -s 0 > log_heap_insert_basic_online_24
 echo "basic heap insert for 2^24"

 
 ./run-experiment -p m:26 c:26 > log_heap_insert_basic_pre_26
 ./run-experiment heap -m 27 -d 26 -i 1 -e 0 -opt 0 -s 0> log_heap_insert_basic_online_26
 echo "basic heap insert for 2^26"

