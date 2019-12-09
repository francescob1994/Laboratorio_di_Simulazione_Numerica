#!/bin/bash

cp Primes Primes_old
sed -i '1d' Primes
./Monte_Carlo_NVT.exe
for i in {1..29}
do
sed -i '1d' Primes
./Monte_Carlo_NVT.exe
done
cp Primes_old Primes
rm Primes_old

#L=(10 20 100 500 1000 5000)
#for i in ${L[@]}
#do
#
#sed -i '6,7d' input.dat
#awk 'NR==6{ print 10000/'$i' }1' input.dat > tmp.dat
#awk 'NR==7{ print '$i' }1' tmp.dat > input.dat
#./Monte_Carlo_NVT.exe
#
#
#done

