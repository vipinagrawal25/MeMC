import numpy as np
import matplotlib.pyplot as plt
import lib as lib
import numpy.linalg as LA
import sys
import os
from collections import defaultdict
from collections import Counter

tot_rings = int(sys.argv[2])
outer_radius = 10
def create_arrays(tot_rings=tot_rings):
    num_fold = 5     # number of connections in disclination
    Nfaces = num_fold*(tot_rings**2)                       # number of triangles  with r rings
    Nedges = int((3*tot_rings**2+tot_rings)*num_fold/2)    # number of edges
    Nverts = int(Nedges-Nfaces+1)                          # number of vertices F-E+V=1 Euler's theorem

    # Vertices
    V=np.zeros((Nverts,3))     # x,y,z
    V[0,:]=[0,0,0]             # the center is at the origin
    count=1
    for r in range(1,tot_rings+1):
        radius = outer_radius * r / tot_rings
        theta=np.linspace(0, 2*np.pi, num_fold*r, endpoint=False)     # divide theta uniformly in num_fold*r parts
        V[count:count+num_fold*r,0] = radius*np.cos(theta)                # x = r cos(theta)
        V[count:count+num_fold*r,1] = radius*np.sin(theta)                # y = r sin(theta)
        V[count:count+num_fold*r,2] = 0                              # z = 0 for a flat case
        count += num_fold*r

    V=V[::-1]
    def distance(i,j): return np.linalg.norm(V[i,:]-V[j,:])

    # Faces
    F=np.zeros((Nfaces,3))                                                     # (v1, v2, v3, +/-1)
    cutoff=outer_radius/tot_rings*np.sqrt(1+(2*np.pi/num_fold)**2)             # lattice length=1
    count=0
    for i in range(Nverts-2):
        for j in range (i+1,Nverts-1):
            if distance(i,j)<cutoff:
                for k in range(j+1,Nverts):
                    if distance(k,i)<cutoff and distance(k,j)<cutoff:
                        # normal=np.cross(V[j,:]-V[i,:],V[k,:]-V[j,:])
                        F[count,:]=[i,j,k]
                        count += 1
    return V,F

def sort_nbrs(R, Np, cmlst, node_nbr):
    zhat = np.array([0.,0.,1.])
    for i in range(Np):
        nbrs=node_nbr[cmlst[i]:cmlst[i+1]]  # neighbours of ith node)
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

def polar(xyz):
    x=xyz[0]
    y=xyz[1]
    z=xyz[2]
    XsqPlusYsq = x**2 + y**2
    return np.arctan2(np.sqrt(XsqPlusYsq),z)

def rotate(vector, nhat, theta):
    '''Rotate a vector about nhat by angle theta using a rotation matrix'''
    nhat = nhat / np.linalg.norm(nhat)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    one_minus_cos = 1 - cos_theta
    nx, ny, nz = nhat
    
    rotation_matrix = np.array([
        [cos_theta + nx**2 * one_minus_cos, nx*ny*one_minus_cos - nz*sin_theta, nx*nz*one_minus_cos + ny*sin_theta],
        [ny*nx*one_minus_cos + nz*sin_theta, cos_theta + ny**2 * one_minus_cos, ny*nz*one_minus_cos - nx*sin_theta],
        [nz*nx*one_minus_cos - ny*sin_theta, nz*ny*one_minus_cos + nx*sin_theta, cos_theta + nz**2 * one_minus_cos]
    ])
    
    rot_vec = np.dot(vector, rotation_matrix.T)
    return rot_vec

def neighbours(Np, simpl, vertices):
    r1=simpl[:,0]
    r2=simpl[:,1]
    lst=np.zeros(Np,dtype=int)
    cumlst=np.zeros(Np+1,dtype=int)
    for i in range(0, Np):
       lst[i]=len(r1[r1==i])/2

    cumlst[1:] = np.cumsum(lst)
    node_neighbour = np.zeros(cumlst[-1],dtype=int)
    for i in range(0, cumlst[-1], 1):
        node_neighbour[i]=r2[2*i]
    ncmlist = np.diff(cumlst)
    # print(node_neighbour)
    return cumlst, new_way_nbrs(Np, vertices[:,0:2], cumlst, node_neighbour, nghst=12)

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
    nsimplices=repeat_unique_pairs(nsimplices)
    nsimplices = np.asarray(sorted(nsimplices, key=lambda x: (x[0], x[1])))
    return nsimplices

def repeat_unique_pairs(arr):
    pairs = [(row[0], row[1]) for row in arr]
    pair_counts = Counter(pairs)
    unique_pairs = [pair for pair, count in pair_counts.items() if count == 1]
    repeated_pairs = np.array([list(pair) + [-1] for pair in unique_pairs])
    return np.vstack([arr, repeated_pairs])

def new_way_nbrs(Np, pts, cmlist, node_nbr, nghst=12):
    new_nbr = np.zeros(nghst*Np, dtype=int)
    new_nbr[:] = -1
    for ip in range(0, Np):
        nbrs = node_nbr[cmlist[ip]:cmlist[ip+1]]
        # print(nbrs)
        num_nbr = -(cmlist[ip]-cmlist[ip+1])
        angles = []
        for i in nbrs:
            dx = pts[i,0]-pts[ip,0]
            dy = pts[i,1]-pts[ip,1]
            angles.append(np.arctan2(dy, dx))
        angles = np.asarray(angles)
        # print(nbrs, angles)
        idx = np.where(angles < 0)
        # angles[idx] = angles[idx] + 2*np.pi
        sort = np.argsort(angles)
        nnbrs = nbrs[sort]
        st_idx = int(ip*nghst); 
        end_idx = int(ip*nghst + num_nbr)
        new_nbr[st_idx:end_idx] = nnbrs[:]
        # print(nbrs,nnbrs)
        # print(new_nbr[st_idx:end_idx])
    return new_nbr

def sort_boundary(Nb, pts, cmlist, node_nbr, nghst=12):
    # ib=3
    # nbrs = node_nbr[nghst*ib:nghst*(ib+1)]
    # print(nbrs)
    # print(pts[ib], pts[nbrs[-1]], pts[nbrs[0:-1]])
    # print(is_all_points_on_one_side(pts[ib], pts[nbrs[-1]], pts[nbrs[0:-1]]))
    for ib in range(Nb):
        num_nbr = -(cmlist[ib]-cmlist[ib+1])
        nbrs = node_nbr[nghst*ib:nghst*ib+num_nbr]
        logic = is_all_points_on_one_side(pts[ib], pts[nbrs[0]], pts[nbrs[1:]]) and is_all_points_on_one_side(pts[ib], pts[nbrs[-1]], pts[nbrs[0:-1]])
        while(logic==False):
            nbrs = np.roll(nbrs, 1)
            logic = is_all_points_on_one_side(pts[ib], pts[nbrs[0]], pts[nbrs[1:]]) and is_all_points_on_one_side(pts[ib], pts[nbrs[-1]], pts[nbrs[0:-1]])
        node_nbr[nghst*ib:nghst*ib+num_nbr] = nbrs
    return node_nbr
    
def is_all_points_on_one_side(vertex1, vertex2, points):
    """
    Check if all points lie on one side of the line joining vertex1 and vertex2.
    """
    line_vec = vertex2 - vertex1
    line_vec /= np.linalg.norm(line_vec)
    
    def side_of_line(point):
        vec = point - vertex1
        return np.cross(line_vec, vec)
    
    sides = np.array([side_of_line(point) for point in points])
    return np.all(sides >= 0) or np.all(sides <= 0)

def reconstruct_neighbor_lists(neighbor_indices, ncmlist):
    """
    Reconstruct a list-of-lists for neighbors.
    For vertex i, its neighbor list is given by a slice of length ncmlist[i]
    from the 1D neighbor_indices array.
    """
    neighbor_lists = []
    start = 0
    for count in ncmlist:
        count = int(count)  # ensure integer
        # Slice out the neighbor list for the current vertex.
        neighbor_lists.append(neighbor_indices[start:start+count].tolist())
        start += count
    return neighbor_lists

def is_boundary_vertex(nbrs, edge_set):
    """
    Check if a vertex is on the boundary. For the given vertex, we examine its neighbor list in consecutive (circular) order. For each consecutive pair (n1, n2), we check if the edge between n1 and n2 exists in the global edge set. If any consecutive pair is missing, we return True (the vertex is on the boundary). Otherwise, we return False (the vertex is interior).
    """
    n = len(nbrs)
    if len(nbrs) < 3:
        # With fewer than two neighbors, we consider it a boundary vertex.
        return True
    # Check every consecutive pair (wrapping from last back to first)
    for i in range(n):
        n1 = int(nbrs[i])
        n2 = int(nbrs[(i + 1) % n])
        edge = (min(n1, n2), max(n1, n2))
        print(edge)
        if edge not in edge_set:
            # Missing the edge connecting two consecutive neighbors => boundary.
            return True
    return False

def make_edge_list(neighbor_indices,ncmlist,nghst=12):
    """
    Construct a unique edge list from a one-dimensional neighbor array.

    Parameters:
        neighbor_indices (array-like): A 1D array containing all the neighbor indices for all vertices.
        ncmlist (array-like): A 1D array where ncmlist[i] is the number of neighbors of vertex i.

    Returns:
        List[Tuple[int, int]]: A sorted list of unique edges as tuples (i, j) with i < j.
    """
    edges = set()  # Use a set to avoid duplicate edges.
    current_index = 0
    # print(neighbor_indices)
    num_vertices = int(len(neighbor_indices)/nghst)
    
    for i in range(num_vertices):
        # Get the neighbor list for vertex i from the flat neighbor_indices array.
        current_index = i*nghst
        for j in neighbor_indices[current_index: current_index + ncmlist[i]]:
            # Create an edge with sorted vertex indices so that (i, j) is same as (j, i)
            edge = (min(i, int(j)), max(i, int(j)))
            edges.add(edge)
    # Return a sorted list of edges for consistency.
    return sorted(list(edges))

vertices,faces=create_arrays()
Np = vertices.shape[0]
print(f"Number of mesh points: {Np}")
print(f"Disc area: {np.pi * 10**2:.3f}")
sort_tri = sort_simplices(faces)
cmlist, node_nbr = neighbours(Np, sort_tri, vertices=vertices)
ncmlist = np.diff(cmlist)
edge_set=make_edge_list(node_nbr,ncmlist)
node_nbr=sort_boundary(5*tot_rings, vertices[:,0:2], cmlist, node_nbr, nghst=12)

lib.write_hdf5(vertices, ncmlist, node_nbr, sys.argv[1])

# bdry_ind = [is_boundary_vertex(node_nbr[12*i:12*i+ncmlist[i]], edge_set) for i in range(Np)]
# print("Boundary vertex :", np.linspace(0,Np-1,Np,dtype=int)[bdry_ind])

# plt.figure(figsize=(6, 6))
# for i, txt in enumerate(ncmlist):
#     plt.annotate(i, (vertices[i, 0], vertices[i, 1]), fontsize=8, ha='right')
# plt.triplot(vertices[:, 0], vertices[:, 1], faces, 'k-', alpha=0.6)
# plt.scatter(vertices[:, 0], vertices[:, 1], color='purple', s=5)
# plt.title("Triangulation with Five-Fold Disclination")
# plt.axis("equal")
# plt.show()

# neighbor_lists = reconstruct_neighbor_lists(node_nbr, ncmlist)
# boundary_vertices = []
# for i in range(len(neighbor_lists)):
#     if is_boundary_vertex(i, neighbor_lists):
#         boundary_vertices.append(i)

# print("Boundary vertex indices:", boundary_vertices)
    
# # Optionally, print their positions.
# print("\nBoundary vertex positions:")
# for i in boundary_vertices:
#     print(f"Vertex {i}: {positions[i]}")


# boundary_vertices = find_boundary_vertices(node_nbr,ncmlist)
# print("Boundary vertex indices:")
# print(sorted(boundary_vertices))



# Plot the triangulation
