import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import lib
import sys
import h5py

def generate_equidistant_points_2d(x_range, y_range, num_points_x, num_points_y):
    x_min, x_max = x_range
    y_min, y_max = y_range

    # Create 1D arrays of x and y coordinates
    x_coords = np.linspace(x_min, x_max, num_points_x)
    y_coords = np.linspace(y_min, y_max, num_points_y)

    # Create a 2D grid of points
    xx, yy = np.meshgrid(x_coords, y_coords)

    # Stack x and y coordinates into a single array
    points = np.column_stack((xx.ravel(), yy.ravel()))

    return points

def generate_vertex_neighbors(faces, num_vertices):
    vertex_neighbors = [set() for _ in range(num_vertices)]
    for face in faces:
        v0, v1, v2 = face
        # Each vertex in a triangle is a neighbor of the other two
        vertex_neighbors[v0].add(v1)
        vertex_neighbors[v0].add(v2)
        vertex_neighbors[v1].add(v0)
        vertex_neighbors[v1].add(v2)
        vertex_neighbors[v2].add(v0)
        vertex_neighbors[v2].add(v1)
    
    # Create the first array: number of neighbors for each vertex
    num_neighbors = np.array([len(neighbors) for neighbors in vertex_neighbors], dtype=np.int32)
    
    # Create the second array: all neighbors concatenated in 1D
    all_neighbors = []
    for neighbors_set in vertex_neighbors:
        all_neighbors.extend(sorted(list(neighbors_set)))
    
    all_neighbors = np.array(all_neighbors, dtype=np.int32)
    
    return num_neighbors, all_neighbors

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

    ncmlist = np.diff(cumlst)

    return cumlst,node_neighbour

def write_hdf5(R, cmlst, node_nbr,  posfile):
    if posfile.split(".")[-1]=="h5":
        pass
    else:
        posfile=posfile+".h5"
    hf = h5py.File(posfile,'w')
    hf.create_dataset('pos',data=R.reshape(-1))
    hf.create_dataset('cumu_list',data=cmlst.astype(np.int32))
    hf.create_dataset('node_nbr',data=node_nbr.astype(np.int32))
    hf.close()

def read_simplices(file):
    simplices = []
    data = np.loadtxt(file, dtype=int)
    for it in range(0,len(data)):
        simplices.append([data[it,0], data[it,1], data[it,2]])
    return simplices

def new_way_nbrs(cmlist, node_nbr, nghst=12):
    new_nbr = np.zeros(nghst*Np, dtype=int)
    new_nbr[:] = -1
    for ip in range(0, Np):
        nbrs = node_nbr[cmlist[ip]:cmlist[ip+1]]
        num_nbr = -(cmlist[ip]-cmlist[ip+1])
        angles = []
        for i in nbrs:
            dx = pts[i,0]-pts[ip,0]
            dy = pts[i,1]-pts[ip,1]
            angles.append(np.arctan2(dy, dx))
        angles = np.asarray(angles)
        idx = np.where(angles < 0)
        angles[idx] = angles[idx] + 2*np.pi
        sort = np.argsort(angles)
        nnbrs = nbrs[sort]
        st_idx = int(ip*nghst); end_idx = int(ip*nghst + num_nbr)
        new_nbr[st_idx:end_idx] = nnbrs[:]
    return new_nbr

# Example usage
if __name__ == "__main__":
    x_range = (-10, 10)  # x-axis range
    y_range = (-10, 10)  # y-axis range
    num_points_x = 16  # Number of points along the x-axis
    num_points_y = 16  # Number of points along the y-axis
    ng=12  # max number of neighbors

    pts = generate_equidistant_points_2d(x_range, y_range, num_points_x, num_points_y)
    Np = len(pts)
    pts_3d = np.zeros(shape=(Np, 3), dtype=float)
    pts_3d[:,0] = pts[:,0]
    pts_3d[:,1] = pts[:,1]
    tri = Delaunay(pts, furthest_site = False)
    faces = tri.simplices
    vertices = pts_3d

    # Use the new function to generate vertex neighbors
    ncmlist, node_nbr = generate_vertex_neighbors(faces, Np)
    cmlist = np.zeros(Np+1, dtype=int)
    cmlist[1:] = np.cumsum(ncmlist)
    new_nbr = new_way_nbrs(cmlist, node_nbr, nghst=ng)
    write_hdf5(vertices,  ncmlist, new_nbr, sys.argv[1])


    # ncmlist = np.diff(cmlist)

    # print("Index | num_neighbors | cmlist")
    # for i in range(Np):
    #     print(f"{i:5d} | {num_neighbors[i]:13d} | {cmlist[i+1] - cmlist[i]:5d}")

    # # exit(1)
    # # Plot the mesh with number of neighbours at each node
    # plt.figure(figsize=(8, 8))
    # plt.triplot(pts[:, 0], pts[:, 1], faces, color='gray', alpha=0.5)
    # plt.scatter(pts[:, 0], pts[:, 1], c=ncmlist, cmap='viridis', s=40, edgecolors='k')
    # for i, (x, y) in enumerate(pts):
    #     plt.text(x, y, str(ncmlist[i]), color='red', fontsize=8, ha='center', va='center')
    # plt.colorbar(label='Number of neighbours')
    # plt.title('Mesh with Number of Neighbours at Each Node')
    # plt.xlabel('X')
    # plt.ylabel('Y')
    # plt.axis('equal')
    # plt.tight_layout()
    # plt.show()

    # nn_nbr = np.zeros(ng*Np, dtype=int)
    # nn_nbr[:] = -1
    # for i in range(0,Np):
    #     cidx = ng*i
    #     k = cidx
    #     for j in range(0,ng-1):
    #         k1 = cidx + j
    #         if(new_nbr[k1+1] != new_nbr[k1]):
    #             nn_nbr[k] = new_nbr[k1];
    #             k = k + 1

    # ncml = np.zeros(Np, dtype = int)
    # ncml[:] = 0
    # for i in range(0,Np):
    #     cidx = ng*i
    #     for j in range(0,ng):
    #         k = cidx + j
    #         if(nn_nbr[k] != -1): ncml[i] = ncml[i] + 1;
