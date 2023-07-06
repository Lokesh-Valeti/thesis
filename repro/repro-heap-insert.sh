  ./run-experiment -p m:16 c:16 > basic_heap_insert_preproc_16
  ./run-experiment -t 64 heap -m 17 -d 16 -i 1 -e 0 -opt 0 -s 0 > basic_heap_insert_online_16
  echo "basic heap insert for 2^16"
 
  ./run-experiment -p m:18 c:18 > basic_heap_insert_preproc_18
  ./run-experiment -t 64 heap -m 19 -d 18 -i 1 -e 0 -opt 0 -s 0 > basic_heap_insert_online_18
  echo "basic heap insert for 2^18"
 
  ./run-experiment -p  m:20 c:20 > basic_heap_insert_preproc_20
  ./run-experiment heap -m 21 -d 20 -i 1 -e 0 -opt 0 -s 0 > basic_heap_insert_online_20
  echo "basic heap insert for 2^20"

 
  ./run-experiment -p m:22 c:22 > basic_heap_insert_basic_preproc_22
  ./run-experiment  heap -m 23 -d 22 -i 1 -e 0 -opt 0 -s 0 > basic_heap_insert_online_22
  echo "basic heap insert for 2^22"

  ./run-experiment -p m:24 c:24 > basic_heap_insert_preproc_24
  ./run-experiment heap -m 25 -d 24 -i 1 -e 0 -opt 0 -s 0 > basic_heap_insert_online_24
  echo "basic heap insert for 2^24"

 
  ./run-experiment -p m:26 c:26 > basic_heap_insert_preproc_26
  ./run-experiment heap -m 27 -d 26 -i 1 -e 0 -opt 0 -s 0> basic_heap_insert_online_26
  echo "basic heap insert for 2^26"



  ./run-experiment -p T0 m:32 r6:1 i4:1 c:5  > opt_heap_insert_preproc_16
  ./run-experiment heap -m 17 -d 16 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_16
  echo "basic heap insert for 2^16"

  ./run-experiment -p T0 m:36 r6:1 i4:1 c:5  > opt_heap_insert_preproc_18
  ./run-experiment heap -m 19 -d 18 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_18
  echo "basic heap insert for 2^18"

  ./run-experiment -p T0 m:40 r6:1 i4:1 c:5  > opt_heap_insert_preproc_20
  ./run-experiment heap -m 20 -d 21 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_20
  echo "basic heap insert for 2^20"

  ./run-experiment -p T0 m:44 r6:1 i4:1 c:5  > opt_heap_insert_preproc_22
  ./run-experiment heap -m 23 -d 22 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_22
   echo "basic heap insert for 2^22"	

  ./run-experiment -p T0 m:48 r6:1 i4:1 c:5  > opt_heap_insert_preproc_24
  ./run-experiment heap -m 25 -d 24 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_24
  echo "basic heap insert for 2^24"

  ./run-experiment -p T0 m:52 r6:1 i4:1 c:5  > opt_heap_insert_preproc_26
  ./run-experiment heap -m 27 -d 26 -i 1 -e 0 -opt 1 -s 0 > opt_heap_insert_online_26
  echo "basic heap insert for 2^26"
