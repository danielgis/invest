import os
from scipy import io
 
irrigation_dir_uri = "/home/mlacayo/maxirrsummarymaps/"
ag_converted_uri = "/home/mlacayo/workspace/ag_converted/Irrigation/"

irrigation_name = "%s_totalmaxirr75.asc"

for base_name in os.listdir(irrigation_dir_uri):
    raster_uri = os.path.join(irrigation_dir_uri, base_name)
    try:
        data = io.loadmat(raster_uri)['totalmaxirr_75areaconstraint']
        
        if data.shape != (4320, 2160):
            raise ValueError
        else:
            crop = base_name[:-13]
            print "Irrigation for %s found." % crop

            asc_uri = os.path.join(ag_converted_uri,
                                   irrigation_name % crop)

            asc_file = open(asc_uri, 'w')
            asc_file.write("ncols        4320\n")
            asc_file.write("nrows        2160\n")
            asc_file.write("xllcorner  -180.0000\n")
            asc_file.write("yllcorner  -90.00000\n")
            asc_file.write("cellsize  8.3333001E-02\n")
            asc_file.write("NODATA_value  -9999.000")

            _, y = data.shape
            for i in range(y):
                asc_file.write("\n"+" ".join([str(v) for v in data[:,i]]).replace("nan", "-9999.00"))

            asc_file.close()
            
    except ValueError:
        print raster_uri, "not readable"
