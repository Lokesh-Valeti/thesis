./run-experiment -p m:17 r17:17 c:17 p:64 > basic_bs_preproc_16
./run-experiment -t 64 bsearch 16 > basic_bs_online_16
echo "[Basic Binary Search for 2^16 done]"

./run-experiment -p m:19 r19:19 c:19 p:64 > basic_bs_preproc_18
./run-experiment -t 64 bsearch 18 > basic_bs_online_18
echo "[Basic Binary Search for 2^18 done]"

./run-experiment -p m:21 r21:21 c:21 p:64 > basic_bs_preproc_20
./run-experiment -t 64 bsearch 20 > basic_bs_online_20
echo "[Basic Binary Search for 2^20 done]"

./run-experiment -p m:23 r23:23 c:23 p:64 > basic_bs_preproc_22
./run-experiment -t 64 bsearch 22 > basic_bs_online_22
echo "[Basic Binary Search for 2^22 done]"

./run-experiment -p m:25 r25:25 c:25 p:64 > basic_bs_preproc_24
./run-experiment -t 64 bsearch 24 > basic_bs_online_24
echo "[Basic Binary Search for 2^24 done]"

./run-experiment -p m:27 r27:27 c:27 p:64 > basic_bs_preproc_26
./run-experiment -t 64 bsearch 26 > basic_bs_online_26
echo "[Basic Binary Search for 2^26 done]"

echo "basic binary search online complete"

./run-experiment -p  i16:1 c:17 p:64 > opt_bs_preproc_16
./run-experiment -t 64 bbsearch 16 > opt_bs_online_16
echo "[Opt Binary Search for 2^16 done]"

./run-experiment -p  i18:1 c:19 p:64 > opt_bs_preproc_18
./run-experiment -t 64 bbsearch 18 > opt_bs_online_18
echo "[Opt Binary Search for 2^18 done]"

./run-experiment -p  i20:1 c:21 p:64 > opt_bs_preproc_20
./run-experiment -t 64 bbsearch 20 > opt_bs_online_20
echo "[Opt Binary Search for 2^20 done]"

./run-experiment -p  i22:1 c:23 p:64 > opt_bs_preproc_22
./run-experiment -t 64 bbsearch 22 > opt_bs_online_22
echo "[Opt Binary Search for 2^22 done]"

./run-experiment -p  i24:1 c:25 p:64 > opt_bs_preproc_24
./run-experiment -t 64 bbsearch 24 > opt_bs_online_24
echo "[Opt Binary Search for 2^24 done]"

./run-experiment -p  i26:1 c:27 p:64 > opt_bs_preproc_26
./run-experiment -t 64 bbsearch 26 > opt_bs_online_26
echo "[Opt Binary Search for 2^26 done]"
~                     
