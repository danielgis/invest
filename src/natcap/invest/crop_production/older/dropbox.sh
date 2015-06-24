#!/bin/bash

SOURCE=~/tmp_source
TMP=~/tmp_destination
DESTINATION=~/ag_converted

rm -rf $SOURCE
rm -rf $TMP
rm -rf $DESTINATION

mkdir $TMP
mkdir $DESTINATION

echo "Setting up folders for OECD income data."
ln -s ~/Dropbox/data\ for\ model/ $SOURCE

echo "Copying compressed source."
cp $SOURCE/OECDIncomeData.nc.gz $TMP

f=$TMP/OECDIncomeData.nc.gz

echo "Expanding compressed source."
gunzip $f

f=$TMP/OECDIncomeData.nc

echo "Converting to TIFF."
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION
rm $f

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

CLIMATE=ClimateBinAnalysis

echo "Setting up folders for climate."
ln -s ~/Dropbox/data\ for\ model/ClimateBinAnalysis/ $SOURCE
mkdir $DESTINATION/$CLIMATE

echo "Copying compressed source."
for f in $SOURCE/*.gz
do
cp $f $TMP
done

echo "Expanding compressed source."
for f in $TMP/*.gz
do
gunzip $f
done

echo "Converting to TIFF."
for f in $TMP/*.nc
do
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$CLIMATE
rm $f
done

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

CROPRANK=CropRank

echo "Setting up folders for croprank."
ln -s ~/Dropbox/data\ for\ model/CropRank $SOURCE
mkdir $DESTINATION/$CROPRANK

echo "Copying compressed source."
for f in $SOURCE/*.gz
do
cp $f $TMP
done

echo "Expanding compressed source."
for f in $TMP/*.gz
do
gunzip $f
done

echo "Converting to TIFF."
for f in $TMP/*.nc
do
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$CROPRANK
rm $f
done

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

FERTILIZER=Fertilizer

echo "Setting up folders for fertilizer."
ln -s ~/Dropbox/data\ for\ model/Fertilizer2000toMarijn $SOURCE
mkdir $DESTINATION/$FERTILIZER

echo "Copying compressed source."
for f in $SOURCE/*.gz
do
cp $f $TMP
done

echo "Expanding compressed source."
for f in $TMP/*.gz
do
gunzip $f
done

echo "Converting to TIFF."
for f in $TMP/*.nc
do
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $f ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$FERTILIZER
rm $f
done

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

MONFREDA=Monfreda

echo "Setting up folders for Monfreda maps."
ln -s ~/Dropbox/data\ for\ model/Monfreda\ maps $SOURCE
mkdir $DESTINATION/$MONFREDA

echo "Expanding compressed source."
for f in $SOURCE/*.zip
do
unzip $f -d $TMP
done

echo "Converting to TIFF."
for f in $TMP/*.asc
do
#fix .asc headers
head -n +6 $f | sed -e 's/^[ \t]*//' > $TMP/tmp.asc
tail -n +7 $f >> $TMP/tmp.asc
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $TMP/tmp.asc ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$MONFREDA
rm $f
done

rm $TMP/tmp.asc

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

NUTRIENT=Nutrient

echo "Setting up folders for nutrient maps."
mkdir $DESTINATION/$MONFREDA/$NUTRIENT

echo "Extracting nutrient maps"
unzip ~/Dropbox/data\ for\ model/Monfreda\ maps/nutrient\ maps.zip -d $TMP

ln -s $TMP/nutrient\ maps $SOURCE

echo "Converting to TIFF."
for f in $SOURCE/*.asc
do
#fix .asc headers
head -n +6 $f | sed -e 's/^[ \t]*//' > $TMP/tmp.asc
tail -n +7 $f >> $TMP/tmp.asc
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $TMP/tmp.asc ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$MONFREDA/$NUTRIENT
rm $f
done

rm $TMP/tmp.asc

rm -rf $SOURCE
rm -rf $TMP

mkdir $TMP

GROUP=11CropGroups

echo "Setting up folders for croprank."
ln -s ~/Dropbox/data\ for\ model/Monfreda\ maps/11CropGroups $SOURCE
mkdir $DESTINATION/$MONFREDA/$GROUP

echo "Expanding compressed source."
for f in $SOURCE/*.zip
do
unzip $f -d $TMP
done

echo "Converting to TIFF."
for f in $TMP/*.asc
do
#fix .asc headers
head -n +6 $f | sed -e 's/^[ \t]*//' > $TMP/tmp.asc
tail -n +7 $f >> $TMP/tmp.asc
gdal_translate -of GTiff -co "COMPRESS=LZW" -a_srs EPSG:4326 $TMP/tmp.asc ${f%.*}.tif
mv ${f%.*}.tif $DESTINATION/$MONFREDA/$GROUP
rm $f
done

rm $TMP/tmp.asc

rm -rf $SOURCE
rm -rf $TMP
