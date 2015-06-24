import os
import re

data_dir_uri = "/home/mlacayo/workspace/data/CropProduction/input"
output_dir_uri = "/home/mlacayo/workspace/guide/source/crop_production"

cbi_mod_yield_dir_uri = os.path.join(data_dir_uri, "CBI_mod_yield")

cbi_mod_yield_pattern = "([A-Z]+[_][a-z]+[_])([a-z]+)"

fertilizer_name = "%s_fertilizer.csv"
statistics_name = "%s_statistics.csv"

unknown_file = []
for base_name in os.listdir(cbi_mod_yield_dir_uri):
    csv_uri = os.path.join(cbi_mod_yield_dir_uri, base_name)
    if os.path.isfile(csv_uri):
        m = re.search(cbi_mod_yield_pattern, base_name)
        if m != None:
            crop = m.group(2)

            fertilizer_uri = os.path.join(output_dir_uri,
                                          fertilizer_name % crop)

            statistics_uri = os.path.join(output_dir_uri,
                                          statistics_name % crop)

            table = [line.split(",") for line in open(csv_uri).readlines()]

            #create fertilizer file
            fertilizer_file = open(fertilizer_uri, 'w')
            for row in table:
                fertilizer_file.write(",".join(row[:13]) + "\n")
            fertilizer_file.close()

            #create staistics table
            statistics_file = open(statistics_uri, 'w')
            for row in table:
                statistics_file.write(",".join(row[:1]+row[21:]).strip().strip(",")+"\n")
            statistics_file.close()

        else:
            unknown_files.append(raster_uri)
