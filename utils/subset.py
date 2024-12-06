import os
import sys
import glob as glob
from visit_utils import *

#files = sorted(glob.glob('/data/vikashp/paris_sing_bub/
# Set a better view
infl = sys.argv[1] #'../main/data/front_0001.silo'
ofl = sys.argv[2]
ResetView()
v = GetView3D()
v.viewNormal = (0.08914879958132636, 0.6611402233492748, 0.744947042817729)
v.viewUp = (-0.03341858151947904, -0.745518469617929, 0.6656466103479887)
v.focus = (0, 0, 0)
v.viewAngle =  30
v.parallelScale = 5.0
v.nearPlane =  -11.5
v.farPlane = 11.5
v.imagePan =  (0, 0)
v.imageZoom =  1.5
v.perspective = 1
v.eyeAngle =  2
SetView3D(v)

ats = AnnotationAttributes()
ats.axes3D.visible = 0
ats.axes3D.triadFlag = 0
ats.axes3D.bboxFlag = 0
ats.userInfoFlag = 0
ats.databaseInfoFlag = 0
ats.legendInfoFlag = 0
ats.axes3D.xAxis.grid = 0
ats.axes3D.yAxis.grid = 0
ats.axes3D.zAxis.grid = 0
ats.axes3D.yAxis.label.font.bold = 0
ats.axes3D.yAxis.label.font.scale = 1
ats.axes3D.setBBoxLocation = 0
# ats.axes3D.bboxLocation = (n1l, n1u, n2l, n2u, n3l, n3u) # figure a way to pass this
SetAnnotationAttributes(ats)

s = SaveWindowAttributes()
#format = PNG, BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
s.format = s.PNG
s.fileName = ofl
s.outputDirectory = "./"
s.width, s.height = 2700,1920
s.screenCapture = 0
SetSaveWindowAttributes(s)

OpenDatabase(infl)
AddPlot("Mesh", "mesh")
AddPlot("Subset", "domains")
sc = SubsetAttributes()
sc.colorType= sc.ColorBySingleColor
#sc.singleColor = (52, 72, 94, 255) # 
sc.singleColor = (149, 165, 166, 255)
SetPlotOptions(sc)
DrawPlots() # Draw the plots in case there is ever an error
SaveWindow()
DeleteAllPlots()
CloseDatabase(infl)
exit()

#( 76, 114, 176, 255)
#(221, 132,  82, 255)
#( 85, 168, 104, 255)
#(196,  78,  82, 255)
#(129, 114, 179, 255)
#(147, 120,  96, 255)
# (231, 76, 60, 255)
# (52, 73, 94, 255)

# Seaborn Paired colors
#(166, 206, 227, 255)
#( 31, 120, 180, 255)
#(178, 223, 138, 255)
#( 51, 160,  44, 255)
#(251, 154, 153, 255)
#(227,  26,  28, 255)
#(253, 191, 111, 255)
#(255, 127,   0, 255)
#(202, 178, 214, 255)
#(106,  61, 154, 255)
#(255, 255, 153, 255)
#(177,  89,  40, 255)
