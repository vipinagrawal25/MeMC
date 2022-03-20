import numpy as np
#
def vtk_points(infile, points, triangles):
    """
    Snippet to write data of a mesh in a specific
    format to be able to load in visit 
    Download visit binaries from
    https://wci.llnl.gov/simulation/computer-codes/visit/executables
    To visualize
    load <output>.vtk in visit
    in Add select subsets>domains
    click draw in gui
    """
    Np = np.shape(points)[0]
    num_triangles = np.shape(triangles)[0]
    ## write the headers of the file
    with open(infile, "w") as f:
        f.write('# vtk DataFile Version 2.0 \n')
        f.write('grid, time 110\n')
        f.write('ASCII \n')
        f.write('DATASET POLYDATA \n')
        f.write('POINTS  '+str(Np)+'  float\n')
        for p in points:
            f.write("%16.8f %16.8f %16.8f\n" %(p[0], p[1], p[2]))
        f.write('POLYGONS  '+str(num_triangles)+'  '
                +str(4*num_triangles) + '\n')
        for it, tri in enumerate(triangles):
            f.write("%d %d %d %d\n" %(3, tri[0], tri[1], tri[2]))
def vtk_points_scalar(infile, points, scalar, name_scalar='points_data'):
    """
    subroutine to dump the scalars evaluated along points
    """
    Np = np.shape(points)[0]
    Nscalar = np.shape(scalar)[0]

    if(Np != Nscalar):
        print ("Error:: the dimension of points and scalar should be same") 
    else:
       ## write the headers of the file
        with open(infile, "a") as f:
            f.write('POINT_DATA  '+str(Np) + '\n')
            f.write('SCALARS '+name_scalar+' float 1 \n')
            f.write('LOOKUP_TABLE default \n')
            for p in scalar:
                f.write("%.16E\n" %(p))
def cells_cells_scalar(infile,
        triangles, 
        scalar, 
        name_scalar='cell_data'):
    """
    subroutine to dump the scalars evaluated along points
    """
    num_triangles = np.shape(triangles)[0]
    Nscalar = np.shape(scalar)[0]

    if(num_triangles != Nscalar):
        print(num_triangles, Nscalar)
        print ("Error:: the number of scalars must match the number of triangles") 
    else:
       ## write the headers of the file
        with open(infile, "a") as f:
            f.write('CELL_DATA  '+str(num_triangles) + '\n')
            f.write('SCALARS '+name_scalar+' float 1 \n')
            f.write('LOOKUP_TABLE default \n')
            for p in scalar:
                f.write("%16.8f\n" %(p))
#

