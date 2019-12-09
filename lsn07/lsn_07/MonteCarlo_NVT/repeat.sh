#!/bin/bash

echo "******************* S I M U L A Z I O N E  *************************"
./Monte_Carlo_NVT.exe > u_p.1
sed -i '1,28d' u_p.1
python correlazioni.py 1

for i in {2..60}
do
echo "******************* S I M U L A Z I O N E  *************************"
sed -i '1d' Primes #cancello la prima riga di primes così la simulazione che segue legge p1 e p2 dalla seconda riga.(seed diverso)
./Monte_Carlo_NVT.exe > u_p.$i
sed -i '1,28d' u_p.$i
python correlazioni.py $i
done

cp Primes_old Primes   #ripristino Primes così come era all'inizio
