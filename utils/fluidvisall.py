import glob
import sys
import os
import numpy as np
################################################################################
def lastfile(dirname,prefix='',zfill=5,suffix=".h5"):
    ''' The function returns the complete path of last file in the 
    #     directory. Here, we use binary search.'''
    lastF=sorted(glob.glob(dirname+prefix+"*"+suffix))[-1]
    low=int(lastF.split("/")[-1].replace(prefix,"").replace(suffix,""))
    return low
################################################################################
dd=sys.argv[1]
nfiles=lastfile(dd,"snap_")
if len(sys.argv)==2:
    jump=1
else:
    jump=int(sys.argv[2])
for i in range(0,nfiles,jump):
	cmd="python utils/fluidvis.py "
	cmd+=dd.replace("/mc_log","/")+"snap_"+str(i).zfill(5)+".h5"
	print(cmd)
	os.system(cmd)
################################################################################