import numpy as np
from numpy import linalg as LA
import quaternion
import os, sys
from scipy.spatial import ConvexHull
import h5py
    

def read_data(filename):
    pos = h5py.File(filename)["pos"][()]
    pos = np.asarray(pos)
    Np = int(len(pos)/2)
    pts_sph = pos.reshape(Np,2)
    pts_cart = np.asarray([[np.sin(theta)*np.cos(phi), 
        np.sin(theta)*np.sin(phi), 
        np.cos(theta)] for theta, phi in pts_sph]
        ) 
    return Np, pts_sph, pts_cart

def triangulate(rr):
    hull = ConvexHull(rr)
    triangles = hull.simplices
    return triangles
##-----------------------------------------------------------------------#

def sort_simplices(cells):
    lsimples = len(cells)
    nsimplices = np.asarray([], dtype=np.int32)
    for scles in cells:
        nscles = np.sort(scles)
        nsimplices = np.hstack([nsimplices, nscles])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[2], nscles[0]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[0], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[0], nscles[2], nscles[1]]])
        nsimplices = np.hstack([nsimplices, [nscles[1], nscles[0], nscles[2]]])
        nsimplices = np.hstack([nsimplices, [nscles[2], nscles[1], nscles[0]]])
    nsimplices = nsimplices.reshape(lsimples*6, 3)
    nsimplices = np.asarray(sorted(nsimplices, key=lambda x: (x[0], x[1])))
    return nsimplices

def neighbours(Np, simpl):
    r1=simpl[:,0]
    r2=simpl[:,1]
    r3=simpl[:,2]
    lst=np.zeros(Np,dtype=int)
    cumlst=np.zeros(Np+1,dtype=int)
    for i in range(0, Np):
       lst[i]=len(r1[r1==i])/2

    cumlst[1:] = np.cumsum(lst)
    node_neighbour = np.zeros(cumlst[-1],dtype=int)
    for i in range(0, cumlst[-1], 1):
        node_neighbour[i]=r2[2*i]
    return cumlst,node_neighbour

def sort_2Dpoints_theta(x,y):
    len_x = len(x)
    len_y = len(y)
    if len_x!=len_y:
        raise Exception("")
    #
    xsort=np.zeros(len_x)
    ysort=np.zeros(len_y)
    #
    theta=np.arctan2(x,y)+np.pi
    indices=np.linspace(0,len_x-1,len_x)
    xyth=np.transpose(np.array([x,y,theta,indices]))
    #
    xysort = np.asarray(sorted(xyth, key=lambda x: (x[2])))
    return xysort[:,3].astype(int),np.array([xysort[:,0],xysort[:,1]])
##-----------------------------------------------------------------------------#

def polar(xyz):
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    XsqPlusYsq = x**2 + y**2
    return np.arctan2(np.sqrt(XsqPlusYsq),z)
##----------------------------------------------------------------------------#

  
def rotate(vector,nhat,theta):
    '''rotate a vector about nhat by angle theta'''
    cos_thby2=np.cos(theta/2)
    sin_thby2=np.sin(theta/2)
    q=np.quaternion(cos_thby2,nhat[0]*sin_thby2,nhat[1]*sin_thby2,nhat[2]*sin_thby2)
    q_inv=np.quaternion(cos_thby2,-nhat[0]*sin_thby2,-nhat[1]*sin_thby2,-nhat[2]*sin_thby2)
    nn=vector.shape[0]
    rot_vec=np.zeros([nn,3])
    for i in range(nn):
        q_vec=np.quaternion(0,vector[i][0],vector[i][1],vector[i][2])
        rot_vec[i]=quater2vec(q*q_vec*q_inv)
    return rot_vec
##----------------------------------------------------------------------------#
def quater2vec(qq,precision=1e-16):
    if qq.w>1e-8:
        print("# ERROR: Quaternion has non-zero scalar value.\n \
               # Can not convert to vector.")
        exit(1)
    return np.array([qq.x,qq.y,qq.z])


#----------------------------------------------------------------------------#
def sort_nbrs(R, Np, cmlst, node_nbr):
    zhat = np.array([0.,0.,1.])
    for i in range(Np):
        nbrs=node_nbr[cmlst[i]:cmlst[i+1]]  # neighbours of ith node
        vector=R[i]
        # I will rotate the coordinate system about this vector
        vhat = np.cross(vector,zhat)       
        vnorm = LA.norm(vhat)
        # If the vector is already lying at z-axis then there is no need to rotate.
        if vnorm>1e-16:
            vhat = vhat/vnorm
            theta = polar(vector)
            # Rotate all the neighbours of a point.
            rotated=rotate(R[nbrs],vhat,theta)
            # Since all the voronoi cells are rotated, sort them in anticlockwise direction
            sorted_indices = sort_2Dpoints_theta(rotated[:,0],rotated[:,1])[0]
            node_nbr[cmlst[i]:cmlst[i+1]]=nbrs[sorted_indices]
    return node_nbr
    #

def write_hdf5(R, cmlst, node_nbr,  cells, file):
    if file.split(".")[-1]=="h5":
        pass
    else:
        file=file+".h5"
    hf = h5py.File(file,'w')
    hf.create_dataset('pos',data=R)
    hf.create_dataset('cumu_list',data=cmlst.astype(np.int32))
    hf.create_dataset('node_nbr',data=node_nbr.astype(np.int32))
    hf.create_dataset('triangles',data=cells.astype(np.int32))
    hf.close()

def write_file(pts_cart, cmlist, node_nbr):
    file = open("../Examples/pts.bin", "wb")
    file.write(pts_cart)
    file.close()

    file = open("../Examples/mesh.bin", "wb")
    file.write(cmlist)
    file.write(node_nbr)
    file.close()


inf = sys.argv[1]
Np, pts_sph, pts_cart = read_data(inf)
triangles = triangulate(pts_cart)
sort_tri = sort_simplices(triangles)
cmlist, node_nbr = neighbours(Np, sort_tri)
node_nbr = sort_nbrs(pts_cart, Np, cmlist, node_nbr)
isDir = os.path.isdir("./conf/") 
if(isDir):
    write_hdf5(pts_cart, cmlist, node_nbr,
             triangles, "./conf/dmemc_conf.h5")
else:
    os.mkdir("conf")
    write_hdf5(pts_cart, cmlist, node_nbr,
             triangles, "./conf/dmemc_conf.h5")

# write_file(pts_cart, cmlist, node_nbr)



