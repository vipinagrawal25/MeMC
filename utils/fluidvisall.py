import glob
import sys
import os
import numpy as np
################################################################################
def lastfile(dirname,prefix='',zfill=5,suffix=".h5"):
    ''' The function returns the complete path of last file in the 
    #     directory. Here, we use binary search.'''
    lastF=sorted(glob.glob(dirname+prefix+"*"+suffix))
    if lastF == []:
        return 0
    else:
        low=int(lastF[-1].split("/")[-1].replace(prefix,"").replace(suffix,""))
    return low
################################################################################
first=0
for dd in sys.argv[1:]:
    # dd=sys.argv[1]
    nfiles=lastfile(dd,"snap_")
    first = lastfile(dd,"snap_",suffix=".vtk")
    first=0
    # if len(sys.argv)==2:
    #     jump=1
    # else:
    #     jump=int(sys.argv[2])
    for i in range(first,nfiles):
        cmd="python utils/fluidvis.py "
        cmd+=dd.replace("/mc_log","/")+"snap_"+str(i).zfill(5)+".h5"
        print(cmd)
        os.system(cmd)
################################################################################