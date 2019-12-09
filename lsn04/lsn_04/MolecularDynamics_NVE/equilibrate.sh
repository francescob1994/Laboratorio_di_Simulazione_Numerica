#!/bin/bash

./clean.sh
sed -i '9d' input.dat
awk 'NR==9{print 0}1' input.dat > tmp.dat
cp tmp.dat input.dat
./MD.exe
sed -i '9d' input.dat
awk 'NR==9{print 1}1' input.dat > tmp.dat
cp tmp.dat input.dat
for i in {1..6}
do
./MD.exe
done

#dopo 10 simulazioni vedo il file ave_temp.dat fino a che temperatura è arrivato. 
#a quel punto faccio clean e faccio un' altra simulazione (con old=1 per non perdere equilibrazione)
#ora su ave_cose.dat ci sono solo le medie dell'ultima simulazione equilibrata. uso questo file per stampare i risultati. 
#dopo 10 esecuzioni vedo che la temperatura oscilla attorno al valore prestabilito in input.dat
