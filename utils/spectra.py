import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.special import legendre
from scipy.special import sph_harm
from scipy.spatial import SphericalVoronoi
import numpy.linalg as la
import math
# from vtk import *
# from vtk.util.numpy_support import vtk_to_numpy
######################## PARAMETERS ###############################
Np=23042
lmax=25
######################### FUNCTIONS ###############################
def vtk_to_data(infile,Np=Np,skiprows=5):
	# Read data from vtk file.
	points=np.loadtxt(infile,skiprows=skiprows,max_rows=Np)
	cells=np.loadtxt(infile,skiprows=skiprows+Np+1,max_rows=2*Np-4)
	return points, cells
# -----------------------------------------------------------------#
def project_onto_sph(points,radius=1):
	''' The function take the points (3 dimensional vector) and 
		projects on an unit sphere. '''
	Np = points.shape[0]
	proj_pnts=np.zeros([Np,3])
	for ip,pos in enumerate(points):
		# x,y,z=pos[0],pos[1],pos[2]
		rad_act=la.norm(pos)
		proj_pnts[ip,0]=pos[0]/rad_act	
		proj_pnts[ip,1]=pos[1]/rad_act
		proj_pnts[ip,2]=pos[2]/rad_act
	return proj_pnts
# -----------------------------------------------------------------#
def cart2sph(points):
	''' Computes theta and phi for the given points.'''
	Np=points.shape[0]
	theta = np.zeros(Np)
	phi = np.zeros(Np)
	for ip,pos in enumerate(points):
		x,y,z=pos[0],pos[1],pos[2]
		#
		xy=np.sqrt(x**2+y**2)
		theta[ip] = np.arctan2(y,x)
		phi[ip] = np.arctan2(xy,z);
	return theta,phi
# -----------------------------------------------------------------#
def height_field_lm(h_theta_phi,theta,phi,area,l,m=0,radius=1):
	m=l
	Np = h_theta_phi.shape[0]
	h_lm=np.zeros(3)
	# theta=list(map(math.degrees,theta))
	# phi=list(map(math.degrees,phi))
	for ip in range(Np):
		h_lm[:] = h_lm[:]+h_theta_phi[ip,:]*sph_harm(m,l,theta[ip],phi[ip]).real*area[ip]\
						  *1/np.sin(theta[ip])*1/radius*1/radius
	return h_lm
########################## SCRIPT ####################################
len_arg = len(sys.argv)
if len_arg>1:
	folder = sys.argv[1]
	filen = sys.argv[2].zfill(5)
else:
	folder="output_fig2_yellow/"
	filen = "02990"
infile = "../"+folder+"/part_"+filen+".vtk"
# infile0 = folder+"/part_00000.vtk"
# -------------------- get the data ------------------------------#
points,cells = vtk_to_data(infile)
# points0,cells0 = vtk_to_data(infile0)
# ----------------------------------------------------------------#
theta,phi = cart2sph(points)
proj_pnts = project_onto_sph(points)
h_theta_phi = points-proj_pnts
sv = SphericalVoronoi(proj_pnts)
area=sv.calculate_areas()
h_lm=np.zeros([lmax,3])
for l in range(0,lmax):
	h_lm[l,:]=height_field_lm(h_theta_phi,theta,phi,area,l)
print(h_lm)
plt.plot(la.norm(h_lm,axis=1),'.-')
plt.show()