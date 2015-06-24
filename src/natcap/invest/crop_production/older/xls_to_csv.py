import os
import xlrd
import unicodecsv

in_uri = '/home/mlacayo/Dropbox/data for model/cropNutrient.xls'
workbook = xlrd.open_workbook(in_uri)
workspace, file_name = os.path.split(in_uri)
base_name, _ = os.path.splitext(file_name)

for sheet_name in workbook.sheet_names():
    out_uri = os.path.join(workspace, "%s_%s.csv" % (base_name, sheet_name))
    out_file = open(out_uri, 'wb')
    csv_writer = unicodecsv.writer(out_file, encoding="utf-8")

    worksheet = workbook.sheet_by_name(sheet_name)
    for i in range(worksheet.nrows):
        row = [cell.value for cell in worksheet.row(i)]        
        csv_writer.writerow(row)

    out_file.close()
