#!/bin/bash

WORKSPACE=~/ag_data/*.gz
OUTPUT=~/ag_data_converted
GZ=~/ag_data_converted/*.gz
NC=~/ag_data_converted/*.nc

rm -rf $OUTPUT
mkdir $OUTPUT

for f in $WORKSPACE
do
cp $f $OUTPUT
done

for f in $GZ
do
gunzip $f
done

for f in $NC
do
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
done
