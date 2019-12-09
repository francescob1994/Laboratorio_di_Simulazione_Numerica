#!/bin/bash
#for i in {0..30}
#do
#for j in {1..8}
#do
#echo "******************* S I M U L A Z I O N E  *************************"
#sed -i '1,2d' input.dat
#awk 'NR==1{print sqrt(5)/2 + ('$i'-15)*sqrt(5./2)*(1-1./sqrt(2))/10 }1' input.dat > tmp.dat
#awk 'NR==1{print '$j'*sqrt(5./2)*(1-1./sqrt(2))/6 }1' tmp.dat > input.dat
#./VMC.exe
#done
#done
for i in {0..20}
do
for j in {0..20}
do
echo "******************* S I M U L A Z I O N E  *************************"
sed -i '1,2d' input.dat
awk 'NR==1{print 0.4 + 0.05*'$i' }1' input.dat > tmp.dat
awk 'NR==1{print 0.1 + '$j'* 0.05 }1' tmp.dat > input.dat
./VMC.exe
done
done
