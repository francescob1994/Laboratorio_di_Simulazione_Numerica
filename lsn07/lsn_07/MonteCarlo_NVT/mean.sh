#!/bin/bash

#cov.mean: <u0u1> <p0p1>
#          <u0u2> <p0p2>
#            .      .
#            .      .

awk '{ print $1, $2 }' cov.1 > tmp.dat
for i in {2..60} #
do
paste tmp.dat cov.$i > tmp2.dat
awk '{ print $1 + $3, $2 + $4 }' tmp2.dat > tmp.dat
done

awk '{print $1/60, $2/60 }' tmp.dat > cov.mean 



#up.mean: <u0> <p0>
#         <u1> <p1>
#          .    .
#          .    .


awk '{ print $1, $2 }' u_p.1 > tmp.dat
for i in {2..60} #
do
paste tmp.dat u_p.$i > tmp2.dat
awk '{ print $1 + $3, $2 + $4 }' tmp2.dat > tmp.dat
done

awk '{print $1/60, $2/60 }' tmp.dat > up.mean
