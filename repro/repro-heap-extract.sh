   ./run-experiment -p m:63 a:16 s:15 r16:90 c:32 p:128 > basic_heap_extract_preproc_16
   echo "preprocessing_heap_16 (basic preproc)"
   
   ./run-experiment heap -m 16 -d 16 -i 0 -e 1 -opt 0 -s 0 > basic_heap_extract_online_16
   echo "preprocessing_heap_16 (basic online)"
 

   echo "preprocessing_heap_16 (opt online)"


   ./run-experiment -p m:71 a:18 s:17 r18:102 c:36 p:128 > basic_heap_extract_preproc_18
   echo "preprocessing_heap_18 (basic preproc)"
   
   ./run-experiment heap -m 18 -d 18 -i 0 -e 1 -opt 0 -s 0 > basic_heap_extract_online_18
   echo "preprocessing_heap_18 (basic online)"

   ./run-experiment -p m:80 a:20 s:19 r20:114 c:40 p:128 > basic_heap_extract_preproc_20
   echo "preprocessing_heap_20 (basic preproc)"
   
   ./run-experiment -t 16 heap -m 20 -d 20 -i 0 -e 1 -opt 0 -s 0 > basic_heap_extract_online_20
   echo "preprocessing_heap_20 (basic online)"


   echo "preprocessing_heap_20 (opt online)"

   ./run-experiment -p m:87 a:22 s:21 r22:126 c:44 p:128 > basic_heap_extract_preproc_22
   echo "preprocessing_heap_22 (basic preproc)"
 
   ./run-experiment -t 16 heap -m 22 -d 22 -i 0 -e 1 -opt 0 -s 0> basic_heap_extract_online_22
   echo "preprocessing_heap_22 (basic online)"
 

   echo "preprocessing_heap_22 (opt online)"

  ./run-experiment -p m:95 a:24 s:23 r24:138 c:48 p:128 > basic_heap_extract_preproc_24
  echo "preprocessing_heap_24 (basic preproc)"
 
  ./run-experiment -t 64 heap -m 24 -d 24 -i 0 -e 1 -opt 0 -s 0> basic_heap_extract_online_24
  echo "preprocessing_heap_24 (basic online)"
 

  echo "preprocessing_heap_24 (opt online)"

  ./run-experiment -p m:103 a:26 s:25 r26:150 c:52 p:128 > basic_heap_extract_preproc_26
  echo "preprocessing_heap_26 (basic preproc)"
 ./run-experiment -t 64 heap -m 26 -d 26 -i 0 -e 1 -opt 0 -s 0 > basic_heap_extract_online_26
  echo "preprocessing_heap_26 (basic online)"
 
  ./run-experiment -p m:103 a:26 s:25 i25.3:1 c:52 p:128 > opt_heap_extract_preproc_26



   ./run-experiment -p  m:63 a:16 s:15  i15.3:1 c:32 p:128 > opt_heap_extract_preproc_16
   echo "preprocessing_heap_16 (opt preproc)"
   
   ./run-experiment heap -m 16 -d 16 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_16


      ./run-experiment -p m:71 a:18 s:17 i17.3:1 c:36 p:128 > opt_heap_extract_preproc_18
   echo "preprocessing_heap_18 (opt preproc)"
   
   ./run-experiment heap -m 18 -d 18 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_18
   echo "preprocessing_heap_18 (opt online)"

      ./run-experiment -p m:80 a:20 s:19 i19.3:1 c:40 p:128 > opt_heap_extract_preproc_20
   echo "preprocessing_heap_20 (opt preproc)"
 
   ./run-experiment heap -m 20 -d 20 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_20
   ./run-experiment -p m:87 a:22 s:21 i21.3:1 c:44 p:128  > opt_heap_extract_preproc_22
   echo "preprocessing_heap_22 (opt preproc)"
 
   ./run-experiment -t 16 heap -m 22 -d 22 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_22
  ./run-experiment -p m:95 a:24 s:23 i23.3:1 c:48 p:128 > opt_heap_extract_preproc_24
  echo "preprocessing_heap_24 (opt preproc)"
 
  ./run-experiment -t 32 heap -m 24 -d 24 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_24
  echo "preprocessing_heap_26 (opt preproc)"
 ./run-experiment -t 64 heap -m 26 -d 26 -i 0 -e 1 -opt 1 -s 0 > opt_heap_extract_online_26
  echo "preprocessing_heap_26 (opt online)"
