#!/bin/bash

SQL=/home/mlacayo/ag_converted/

BASE=/home/mlacayo/ag_converted

CLIMATE=ClimateBinAnalysis
CROPRANK=CropRank
FERTILIZER=Fertilizer

MONFREDA=Monfreda

NUTRIENT=Nutrient
GROUP=11CropGroups

for p in $BASE/$CLIMATE $BASE/$CROPRANK $BASE/$FERTILIZER $BASE/$MONFREDA $BASE/$MONFREDA/$NUTRIENT $BASE/$MONFREDA/$GROUP
do
/usr/local/bin/raster2pgsql -s 4326 -d -I -C -M -R $(/usr/bin/find $p -name "*.tif" | /usr/bin/xargs) > $SQL/${p##*/}.sql
done
