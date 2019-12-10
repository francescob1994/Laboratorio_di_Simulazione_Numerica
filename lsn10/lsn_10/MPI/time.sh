#!/bin/bash
rm mean.dat
for cores in {1..4}
do
sed -i '9d' input.dat
awk 'NR==9{ print '$cores'}1' input.dat > tmp.dat
cp tmp.dat input.dat
for i in {1..50}
do
mpiexec -np $cores ./nTSP.exe > output.dat
awk 'NR==13{ print $1 }' output.dat >> time.dat
done
awk '{ av += $1; av2+= $1*$1 }END{ print '$cores', av/NR, sqrt((av2/NR-(av/NR)*(av/NR))/NR) }' time.dat >> mean.dat
echo "ncores = '$cores'"
echo "time"
cat time.dat
echo "mean"
cat mean.dat
rm -rf time.dat
done
rm tmp.dat
rm -rf output.dat

