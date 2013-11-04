#!/bin/bash
imax=$1
field=$2
if [ "x$imax" = "x" ];then
    echo "Please provide the maximum snapshot number"
    exit 1
fi

if [ "x$field" = "x" ];then
    field=0
fi

rm -rf animation/*.png
for iobs in $(seq 1 1 $imax)
do
    python comet-analysis.py Observation iobs=$iobs,field=$field
done
echo "Animating..."
convert -delay 100 animation/picture-*.png -loop 0 animation/disintegration.gif
