import numpy as np
import glob as glob
import lib as lib
import sys




def make_triangles(start, nbrs, triangles):
    for i in range(start, start+len(nbrs)-1):
        triangles.append([start, i+1, i+2])
    triangles.append([start,start+1,start+len(nbrs)])
    return start + len(nbrs)+1

def triangulate_solids(all_nbrs):
    tri_all = []
    id_n = 0
    for nbr_s in all_nbrs:
        nidn = make_triangles(id_n, nbr_s, tri_all)
        id_n = nidn
    return tri_all



dirn = sys.argv[1]
nghst = 12
dtmp = lib.readHdf5(dirn+"/snap_00000.h5", "pos")
Np = int(len(dtmp)/3)
pts_3d = dtmp.reshape(Np,3)
nframe = int(4*np.sqrt(Np))
solid_file = "/stick_pos_id.h5"

# solids = lib.readHdf5(dirn+solid_file, "stick_id")
# solid_id = np.where(solids == True)[0]
# solid_id = np.asarray([128, 300, 500, 328, 512])
all_pts =  np.full(Np, True, dtype=bool)
all_pts[0:nframe] = False 
all_id = np.where(all_pts == True)[0]

Nnbr = lib.readHdf5(dirn+'/snap_00000.h5', "node_nbr")
Cmlst = lib.readHdf5(dirn+'/snap_00000.h5', "cumu_list")
# num_nbr_solid = Cmlst[solid_id]

# nbr_solid = []
# for k, id_ in enumerate(solid_id):
#     nnbr_id = num_nbr_solid[k]
#     st_idx = nghst*id_
#     nbr_solid.append(Nnbr[st_idx:st_idx+nnbr_id])

# tri_sld = triangulate_solids(nbr_solid)

files = sorted(glob.glob(dirn+'snap_00*.h5'))
#fse = sorted(glob.glob(dirn+'stretch_ener_00*.h5'))
jump = int(sys.argv[2])
ix = 550
for i, f in enumerate(files[1::jump]):
    dtmp = lib.readHdf5(f, "pos")
    Np = int(len(dtmp)/3)
    pts_3d = dtmp.reshape(Np,3)
    Nnbr = lib.readHdf5(f, "node_nbr")
    Cmlst = lib.readHdf5(f, "cumu_list")
    num_nbr_all = Cmlst[all_id]
    nbr_all = []
    for k, id_ in enumerate(all_id):
       nnbr_id = num_nbr_all[k]
       st_idx = 12*id_
       nbr_all.append(Nnbr[st_idx:st_idx+nnbr_id])
    tri_all = triangulate_solids(nbr_all)

    # outf = "%sviz_solid_%04d.vtk" %(dirn,i)
    # lib.vtk_solid(outf, pts_3d, solid_id, nbr_solid, tri_sld)

#     outf = "%ssolidpts_%04d.3D" %(dirn,i)
#     lib.vtkScatter(outf, pts_3d[solid_id])
    outf = "%sviz_all_%04d.vtk" %(dirn,i)
    lib.vtk_solid(outf, pts_3d, all_id, nbr_all, tri_all)

#    scalar = lib.readHdf5(fse[i], 'E')
    # outf = "%straj_%04d.3D" %(dirn,i)
    # f = open(outf, 'w')
    # f.write("#x y z v\n" )
    # f.write("%0.5f %0.5f %0.5f 0.1\n" %(pts_3d[ix,0], pts_3d[ix,1], pts_3d[ix,2]))
    # f.write("%0.5f %0.5f %0.5f 0.1\n" %(pts_3d[ix,0], pts_3d[ix,1], pts_3d[ix,2]))
    # f.close()

    # outf = "%strajr_%04d.3D" %(dirn,i)
    # lib.vtkScatter(outf, np.asarray([pts_3d[0,ix], pts_3d[1,ix], pts_3d[2,ix]]))

#    lib.vtk_points_scalar(outf, pts_3d, scalar, name_scalar='se')
