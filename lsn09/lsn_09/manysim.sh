#!/bin/bash
rm tau.dat

for i in {1..100}
do
./clean.sh
./TSP.exe
python GA.py
done

awk '{ av_min+=$1; av_ave+=$2 }END{print "average tau_min = ", av_min/NR, "\naverage tau_ave = ", av_ave/NR}' tau.dat

