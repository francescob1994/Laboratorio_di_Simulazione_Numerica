#!/bin/bash
rm tau.dat

for i in {1..100}
do
./clean.sh
./TSP.exe
python GA.py
done

awk '{ av_min+=$1; av_ave+=$2 }END{print "average tau_min = ", av_min/100, "\naverage tau_ave = ", av_ave/100}' tau.dat

