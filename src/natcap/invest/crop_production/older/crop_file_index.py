import re
import os
import random

rasters = {}

data_dir_uri = "/home/mlacayo/workspace/data/CropProduction/input"

fertilizer_dir_uri = os.path.join(data_dir_uri, "Fertilizer")
climate_dir_uri = os.path.join(data_dir_uri, "ClimateBinAnalysis")
production_dir_uri = os.path.join(data_dir_uri, "Monfreda")
cbi_dir_uri = os.path.join(data_dir_uri, "CBI")
cbi_yield_dir_uri = os.path.join(data_dir_uri, "CBI_yield")
income_climate_dir_uri = os.path.join(data_dir_uri, "income_climate")
yield_mod_dir_uri = os.path.join(data_dir_uri, "yield_mod")

fertilizer_pattern = "([a-z]+)([0-9A-Z]+)([a-z]+)"
climate_pattern = "([A-Za-z]+[_])([a-z]+)([_])([A-Za-z]+[_][0-9])([_][A-Za-z]+[_][0-9]+[_][A-Za-z]+[_][0-9]+[x][0-9]+[_][a-z]+[_])(BinMatrix)"
production_pattern = "([a-z]+)([_])([a-z]+)"
cbi_pattern = "([A-Z]+)([_][a-z]+[_])([a-z]+)"
cbi_yield_pattern = "([A-Z]+[_][a-z]+)([_])([a-z]+)"
income_climate_pattern = "([A-Z]+[_][a-z]+[_])([a-z]+)"
yield_mod_pattern = "([a-z]+)"

file_index_uri = "/home/mlacayo/workspace/data/CropProduction/input/file_index.csv"


unknown_files = []
column_header = set()
for base_name in os.listdir(fertilizer_dir_uri):
    raster_uri = os.path.join(fertilizer_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(fertilizer_pattern, base_name)
        if m is not None:
            crop = m.group(1)
            column = "_".join([m.group(2), m.group(3)])
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

bin_matrix_columns = set()
for base_name in os.listdir(climate_dir_uri):
    raster_uri = os.path.join(climate_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(climate_pattern, base_name)
        if m is not None:
            crop = m.group(2)
            column = "_".join([m.group(4), m.group(6)])
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
            bin_matrix_columns.add(column)
        else:
            unknown_files.append(raster_uri)

for base_name in os.listdir(production_dir_uri):
    raster_uri = os.path.join(production_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(production_pattern, base_name)
        if m is not None:
            crop = m.group(1)
            column = m.group(3)
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

for base_name in os.listdir(cbi_dir_uri):
    raster_uri = os.path.join(cbi_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(cbi_pattern, base_name)
        if m is not None:
            crop = m.group(3)
            column = m.group(1)
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

for base_name in os.listdir(cbi_yield_dir_uri):
    raster_uri = os.path.join(cbi_yield_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(cbi_yield_pattern, base_name)
        if m is not None:
            crop = m.group(3)
            column = m.group(1)
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

for base_name in os.listdir(income_climate_dir_uri):
    raster_uri = os.path.join(income_climate_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(income_climate_pattern, base_name)
        if m is not None:
            crop = m.group(2)
            column = "Income_Climate"
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

for base_name in os.listdir(yield_mod_dir_uri):
    raster_uri = os.path.join(yield_mod_dir_uri, base_name)
    if os.path.isfile(raster_uri):
        m = re.search(yield_mod_pattern, base_name)
        if m is not None:
            crop = m.group(1)
            column = "Yield_mod"
            column_header.add(column)

            if not (crop in rasters):
                rasters[crop] = {}

            rasters[crop][column] = raster_uri[len(data_dir_uri)+1:]
        else:
            unknown_files.append(raster_uri)

column_header = list(column_header)
column_header.sort()

crops = list(rasters.keys())
crops.sort()

file_index = open(file_index_uri, 'w')
file_index.write(",".join(["Id", "Crop"]+column_header))

for i, crop in enumerate(crops):
    row = [str(i+1), crop]
    for column in column_header:
        try:
            row.append(rasters[crop][column])
        except KeyError:
            row.append("")
    file_index.write("\n"+",".join(row))

file_index.close()


#generate percentile and modeled tabs
bin_matrix_columns = list(bin_matrix_columns)
bin_matrix_crops = []
for i, crop in enumerate(crops):
    raster = False
    for column in bin_matrix_columns:
        try:
            if rasters[crop][column] != "":
                if raster is True:
                    print i + 1, " ", crop, " has multiple BinMatrix rasters"
                    break
                else:
                    raster = True
                    bin_matrix_crops.append(crop)
        except KeyError:
            pass
    if raster is False:
        print i + 1, " ", crop, " has no BinMatrix rasters"

iui_percentile_crop = """
{
        "id": "crop_%s_percentile",
        "args_id": "crop_%s_percentile",
        "type": "checkbox",
        "label":"%s",
        "defaultValue": false,
        "required":false
},"""

for crop in  bin_matrix_crops:
    print iui_percentile_crop % (crop, crop, crop.title()),


##print column_header
##print crops
##print unknown_files
##
##for crop in crops:
##    missing = list(set(column_header).difference(set(rasters[crop].keys())))
##    if len(missing) > 0:
##        missing.sort()
##        print "\nCrop %s missing column(s): %s" % (crop, ", ".join(missing))
