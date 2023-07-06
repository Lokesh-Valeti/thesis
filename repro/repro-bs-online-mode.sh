./run-experiment -o -t 64 bsearch 16 > basic_bs_complete_16
echo "[Basic Binary Search for 2^16 done]"
./run-experiment -o -t 64 bsearch 18 > basic_bs_complete_18
echo "[Basic Binary Search for 2^18 done]"
./run-experiment -o -t 64 bsearch 20 > basic_bs_complete_20
echo "[Basic Binary Search for 2^20 done]"
./run-experiment -o -t 64 bsearch 22 > basic_bs_complete_22
echo "[Basic Binary Search for 2^22 done]"
./run-experiment -o -t 64 bsearch 24 > basic_bs_complete_24
echo "[Basic Binary Search for 2^24 done]"
./run-experiment -o -t 64 bsearch 26 > basic_bs_complete_26
echo "[Basic Binary Search for 2^26 done]"

echo "basic binary search online complete"

./run-experiment -o -t 64 bbsearch 16 > opt_bs_complete_16
echo "[Opt Binary Search for 2^16 done]"
./run-experiment -o -t 64 bbsearch 18 > opt_bs_complete_18
echo "[Opt Binary Search for 2^18 done]"
./run-experiment -o -t 64 bbsearch 20 > opt_bs_complete_20
echo "[Opt Binary Search for 2^20 done]"
./run-experiment -o -t 64 bbsearch 22 > opt_bs_complete_22
echo "[Opt Binary Search for 2^22 done]"
./run-experiment -o -t 64 bbsearch 24 > opt_bs_complete_24
echo "[Opt Binary Search for 2^24 done]"
./run-experiment -o -t 64 bbsearch 26 > opt_bs_complete_26
echo "[Opt Binary Search for 2^26 done]"
~                     
