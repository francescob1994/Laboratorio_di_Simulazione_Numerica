#!/bin/bash
for i in {0..30}
do
echo "******************* S I M U L A Z I O N E  $((i+1)) *************************"
sed -i '1d' input.dat
awk 'NR==1{print 0.5+'$i'*(2-0.5)/30}1' input.dat > tmp.dat
cp tmp.dat input.dat
./Monte_Carlo_ISING_1D.exe
done
