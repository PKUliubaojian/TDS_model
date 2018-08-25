# MakeXYLayer.py
# Description: Creates an XY layer and exports it to a layer file

# import system modules 
import arcpy
import os
# Set environment settings
#这里分别是输入和输出文件夹
arcpy.env.workspace = "D:\data\GNSS\CYGNSS"
output="D:\data\GNSS\CYGNSS2"
ws=arcpy.env.workspace
try:
    files=os.listdir(ws)
    #没有做额外的筛选，只支持.csv文件，不太建议文件夹下有别的文件类型
    for name in files:
    # Set the local variables
        in_Table = name
        #经纬度字段要和csv保持一致
        x_coords = "Lon"
        y_coords = "Lat"
        temp=name.strip(".csv")
        out_Layer = temp+"layer"
        saved_Layer = temp+".shp"
        # Set the spatial reference
        spRef = r"Coordinate Systems\Geographic Coordinate Systems\World\WGS 1984.prj"
        # Make the XY event layer...
        arcpy.MakeXYEventLayer_management(in_Table, x_coords, y_coords, out_Layer,spRef)
        # Print the total rows
        print(arcpy.GetCount_management(out_Layer))
        
        # Save to a layer file
        arcpy.FeatureClassToFeatureClass_conversion(out_Layer,output, saved_Layer)
        
        Input_Table = output+'\\'+saved_Layer
        # Process: Add Field
        arcpy.AddField_management(Input_Table, "SNR2", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")
        # Process: Calculate Field
        arcpy.CalculateField_management(Input_Table, "SNR2", "CDbl([SNR])", "VB", "")


except Exception as err:
    print(err.args[0])

