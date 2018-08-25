import arcgisscripting, sys, string, os
gp = arcgisscripting.create()
try:
    pathin="D:\data\GNSS\CYGNSS"
    pathout="D:\data\GNSS\CYGNSS2"
    files=os.listdir(pathin)
    for name in files:
        in_Table = pathin+"/"+name
        in_x = "x"
        in_y = "y"
        temp=name.strip(".csv")
        out_Layer =temp.replace("-","_")
        spref = r"Coordinate Systems\Geographic Coordinate Systems\World\WGS 1984.prj"
        gp.MakeXYEventLayer(in_Table, in_x, in_y, out_Layer, spref)
        pDSC = gp.describe(out_Layer)
        print gp.getcount(out_Layer)
        gp.FeatureClassToShapefile(out_Layer,pathout)
except:
    print gp.GetMessages()