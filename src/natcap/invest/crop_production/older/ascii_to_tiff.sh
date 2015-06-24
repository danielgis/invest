#!/bin/bash

WORKSPACE=~/ag_data/*.zip
OUTPUT=~/ag_data_converted
ASCII=~/ag_data_converted/*.asc

rm -rf $OUTPUT
mkdir $OUTPUT

for f in $WORKSPACE
do
unzip $f -d $OUTPUT
done

for f in $ASCII
do
#fix .asc headers
head -n +6 $f | sed -e 's/^[ \t]*//' > $OUTPUT/tmp.asc
tail -n +7 $f >> $OUTPUT/tmp.asc
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $OUTPUT/tmp.asc ${f%.*}.tif
done
