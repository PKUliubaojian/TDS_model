# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# str2num.py
# Created on: 2018-05-22 17:39:06.00000
#   (generated by ArcGIS/ModelBuilder)
# Usage: str2num <Input_Table> 
# Description: 
# ---------------------------------------------------------------------------

# Import arcpy module
import arcpy

# Script arguments
Input_Table = arcpy.GetParameterAsText(0)

# Local variables:
Output_Feature_Class = Input_Table
Output_Feature_Class__2_ = Output_Feature_Class

# Process: Add Field
arcpy.AddField_management(Input_Table, "SNR2", "DOUBLE", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", "")

# Process: Calculate Field
arcpy.CalculateField_management(Output_Feature_Class, "SNR2", "CDbl([SNR])", "VB", "")

