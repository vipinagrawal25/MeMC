import numpy as np
def write_param(fname="../para_file.in",paramdict=None):
	if paramdict is None:
		paramdict={"N":5120,"coef_bending":2.5,"Y":25,"coef_vol_expansion":1e6,
					"sp_curv":2,"pressure":0,"radius":1.0,
					"pos_bot_wall":-1.05,"sigma":0.17,"epsilon":4,"th_cr":0.53,
					"algo":"mpolis","Dfac":4,"kBT":1,"is_restart":1,
					"mc_total_iters":10000,"mc_dump_iter":10,
					"type":"random","minA":0.95,"maxA":1.05,
					"tip_radius":0.2,"tip_pos_z":1.05,"afm_sigma":0.17,"afm_epsilon":4,
					"spr_compute":0,"nPole_eq_z":-1,"sPole_eq_z":1}
	with open(fname, "w") as file:
		file.write("## Membrane parameters\n")
		file.write("N coef_bending Y coef_vol_expansion sp_curv pressure radius\n")
		file.write("%d %3.2f %3.2f %3.2f %3.2f %3.2f %3.2f\n" %(paramdict['N'],paramdict['coef_bending'],
					paramdict['Y'],paramdict['coef_vol_expansion'],paramdict['sp_curv'],
					paramdict['pressure'],paramdict['radius']))
		file.write("## Sticking parameters\n")
		file.write("pos_bot_wall sigma epsilon th_cr\n")
		file.write("%3.2f %3.2f %3.2f %3.2f\n" %(paramdict['pos_bot_wall'],
					paramdict['sigma'],paramdict['epsilon'],paramdict['th_cr']))
		file.write("## Montecarlo parameters\n")
		file.write("algo Dfac kBT is_restart mc_total_iters mc_dump_iter\n")
		file.write("%s %d %3.2f %d %d %d\n" %(paramdict['algo'],
					paramdict['Dfac'],paramdict['kBT'],paramdict['is_restart'],
					paramdict['mc_total_iters'],paramdict['mc_dump_iter']))
		#
		file.write("## Activity parameters\n")
		file.write("type minA maxA\n")
		file.write("%s %3.2f %3.2f\n" %(paramdict['type'],
					paramdict['minA'],
					paramdict['maxA']))
		#
		file.write("## Afm Tip parameters\n")
		# file.write("N extent_l extent_r extent_t extent_b\n")
		# file.write("%d %1.1f %1.1f %1.1f %1.1f\n" %(paramdict['afm_N'],paramdict['extent_l'],
					# paramdict['extent_r'],paramdict['extent_t'],paramdict['extent_b']))
		file.write("tip_radius tip_pos_z afm_sigma afm_epsilon\n")
		file.write("%2.2f %2.2f %2.2f %2.2f\n" %(paramdict['tip_radius'],paramdict['tip_pos_z'],
					paramdict['afm_sigma'],paramdict['afm_epsilon']))
		file.write("## Spring parameters\n")
		file.write("spr_compute nPole_eq_z sPole_eq_z\n")
		file.write("%d %2.3f %2.3f\n" %(paramdict['spr_compute'],paramdict['nPole_eq_z'],
					paramdict['sPole_eq_z']))
#
def read_param(fname='../para_file.in'):
	## Membrane parameters
	N,coef_bending,Y,coef_vol_expansion,sp_curv,pressure,radius = np.loadtxt(fname,skiprows=2,max_rows=1)
	N=int(N)
	pos_bot_wall,sigma,epsilon,th_cr = np.loadtxt(fname,skiprows=5,max_rows=1)
	## Montecarlo parameters
	algo,Dfac,kBT,is_restart,mc_total_iters,mc_dump_iter=np.loadtxt(fname,skiprows=8,max_rows=1,
																	dtype=str)
	Dfac=int(Dfac)
	is_restart=int(is_restart)
	mc_total_iters=int(mc_total_iters)
	mc_dump_iter=int(mc_dump_iter)
	kBT=float(kBT)
	## Afm Tip parameters
	# afm_N,extent_l,extent_r,extent_t,extent_b = np.loadtxt(fname,skiprows=10,max_rows=1)
	# afm_N=int(afm_N)
	act_type,minA,maxA = np.loadtxt(fname,skiprows=11,max_rows=1,dtype=str)
	minA,maxA=float(minA),float(maxA)
	tip_radius,tip_pos_z,afm_sigma,afm_epsilon = np.loadtxt(fname,skiprows=14,max_rows=1)
	spr_compute,nPole_eq_z,sPole_eq_z = np.loadtxt(fname,skiprows=17,max_rows=1)
	paramdict={"N":N,"coef_bending":coef_bending,"Y":Y,
				"coef_vol_expansion":coef_vol_expansion,"sp_curv":sp_curv,
				"pressure":pressure,"radius":radius,
				"pos_bot_wall":pos_bot_wall,"sigma":sigma,"epsilon":epsilon,"th_cr":th_cr,
				"algo":algo,"Dfac":Dfac,"kBT":kBT,"is_restart":is_restart,"mc_total_iters":mc_total_iters,
				"mc_dump_iter":mc_dump_iter,
				"type":act_type,"minA":minA,"maxA":maxA,
				"tip_radius":tip_radius,"tip_pos_z":tip_pos_z,"afm_sigma":afm_sigma,
				"afm_epsilon":afm_epsilon,
				"spr_compute":spr_compute,"nPole_eq_z":nPole_eq_z,"sPole_eq_z":sPole_eq_z}
	return paramdict
#
def change_param(finname="para_file.in",foutname=None,**kwargs):
	''' The function change one parameter from the input file and overwrites the new parameter file.'''
	if foutname is None:
		foutname=finname
	# print(finname)
	paramdict=read_param(fname=finname)
	for key,value in kwargs.items():
		paramdict[key]=value
	write_param(fname=foutname,paramdict=paramdict)
#