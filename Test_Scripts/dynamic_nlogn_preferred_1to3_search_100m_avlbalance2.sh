#! /bin/bash
# Create a.out by compiling with -D NLOGN -D ENABLE_PREFERRED_TEST -D ENABLE_1TO3 -D AVL_BALANCE -D HEIGHT_DIFF=2
read -s -p "Enter Password for sudo: " sudoPW
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 62746 -b -g -j -m 100 -i 100 > ./test/dynamic/63k/i7.14.63k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 118650 -b -g -j -m 100 -i 100 > ./test/dynamic/119k/i7.14.119k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 172455 -b -g -j -m 100 -i 100 > ./test/dynamic/172k/i7.14.172k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 224979 -b -g -j -m 100 -i 100 > ./test/dynamic/225k/i7.14.225k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 276589 -b -g -j -m 100 -i 100 > ./test/dynamic/277k/i7.14.277k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 327491 -b -g -j -m 100 -i 100 > ./test/dynamic/327k/i7.14.327k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 377820 -b -g -j -m 100 -i 100 > ./test/dynamic/378k/i7.14.378k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 427667 -b -g -j -m 100 -i 100 > ./test/dynamic/428k/i7.14.428k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 477101 -b -g -j -m 100 -i 100 > ./test/dynamic/477k/i7.14.477k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 526172 -b -g -j -m 100 -i 100 > ./test/dynamic/526k/i7.14.526k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 1003201 -b -g -j -m 100 -i 100 > ./test/dynamic/1003k/i7.14.1003k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 1464689 -b -g -j -m 100 -i 100 > ./test/dynamic/1465k/i7.14.1465k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 1916615 -b -g -j -m 100 -i 100 > ./test/dynamic/1917k/i7.14.1917k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 2361679 -b -g -j -m 100 -i 100 > ./test/dynamic/2362k/i7.14.2362k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 2801418 -b -g -j -m 100 -i 100 > ./test/dynamic/2801k/i7.14.2801k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 3236823 -b -g -j -m 100 -i 100 > ./test/dynamic/3237k/i7.14.3237k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 3668582 -b -g -j -m 100 -i 100 > ./test/dynamic/3669k/i7.14.3669k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 4097202 -b -g -j -m 100 -i 100 > ./test/dynamic/4097k/i7.14.4097k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
echo $sudoPW | sudo -S time taskset -c 0 ./a.out -d 3 -n 4523071 -b -g -j -m 100 -i 100 > ./test/dynamic/4523k/i7.14.4523k.3d.nlogn.preferred.1to3.search.100m.avlbalance2.1t.100i.txt
