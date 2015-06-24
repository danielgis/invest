#!/bin/bash
DIR=~/workspace/ag_converted/Irrigation

for f in $DIR/*.asc
do
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
done
