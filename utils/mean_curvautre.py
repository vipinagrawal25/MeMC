import numpy as np
import sys
import h5py
import glob
import os
import vtk
from vtk.util import numpy_support
from multiprocessing import Pool, cpu_count, Manager, Queue
from datetime import datetime
import argparse

class Mesh:
    def __init__(self, points, faces):
        self.points = points
        self.faces = faces

def read_vtk_file(filename):
    """
    Read a VTK file and return a Mesh object containing points and faces
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    
    polydata = reader.GetOutput()
    
    # Get points
    points_vtk = polydata.GetPoints().GetData()
    points = numpy_support.vtk_to_numpy(points_vtk)
    
    # Get faces
    cells = polydata.GetPolys()
    n_cells = cells.GetNumberOfCells()
    faces = np.zeros((n_cells, 4), dtype=np.int64)
    
    cells.InitTraversal()
    for i in range(n_cells):
        idList = vtk.vtkIdList()
        cells.GetNextCell(idList)
        faces[i,0] = idList.GetNumberOfIds()  # Should be 3 for triangles
        for j in range(faces[i,0]):
            faces[i,j+1] = idList.GetId(j)
    
    return Mesh(points, faces)

def compute_mean_curvature(mesh):
    """
    Compute the mean curvature at each vertex of a triangular mesh.

    The method uses the cotangent formula to approximate the Laplace–Beltrami operator.
    For each vertex i, the Laplace operator is approximated by:
    
        Δx_i = (1 / (2A_i)) * Σ_j (w_ij * (x_j - x_i))
    
    where:
      - A_i is the vertex area (computed as 1/3 of the area of all incident triangles),
      - w_ij = cot(α_ij) + cot(β_ij) are the weights computed from the angles opposite the edge (i,j).
    
    The mean curvature H at vertex i is then defined as:
    
        H(i) = 0.5 * || Δx_i ||
    
    Parameters:
        mesh (pyvista.PolyData): A triangular mesh loaded from a VTK file.
    
    Returns:
        numpy.ndarray: Array of mean curvature values, one per vertex.
    """
    points = mesh.points
    # Assumes that the mesh is triangular.
    # PyVista stores faces in a flat array where each face is prefixed by the number of points.
    faces = mesh.faces.reshape(-1, 4)[:, 1:4]
    N = points.shape[0]

    # Allocate arrays for vertex areas and Laplace contributions.
    vertex_area = np.zeros(N)
    laplace = np.zeros((N, 3))
    
    # Dictionary to accumulate cotan weights for each (undirected) edge.
    edge_weights = {}

    # Helper function to compute the cotangent of the angle between two vectors.
    def cotan(u, v):
        # Avoid division by zero
        cross_norm = np.linalg.norm(np.cross(u, v))
        if cross_norm < 1e-8:
            return 0.0
        return np.dot(u, v) / cross_norm

    # Loop over each triangle face.
    for tri in faces:
        i, j, k = tri
        v_i, v_j, v_k = points[i], points[j], points[k]
        
        # Compute edge vectors.
        e_ij = v_j - v_i
        e_ik = v_k - v_i
        e_jk = v_k - v_j
        
        # Triangle area (using half the norm of the cross product).
        tri_area = 0.5 * np.linalg.norm(np.cross(e_ij, e_ik))
        
        # Distribute 1/3 of the triangle area to each vertex.
        vertex_area[i] += tri_area / 3.0
        vertex_area[j] += tri_area / 3.0
        vertex_area[k] += tri_area / 3.0
        
        # In a triangle, the angles are as follows:
        # - Angle at vertex i is between e_ij and e_ik, opposite edge (j, k)
        # - Angle at vertex j is between -e_ij and e_jk, opposite edge (i, k)
        # - Angle at vertex k is between -e_ik and -e_jk, opposite edge (i, j)
        cot_angle_at_i = cotan(e_ij, e_ik)
        cot_angle_at_j = cotan(-e_ij, e_jk)
        cot_angle_at_k = cotan(-e_ik, -e_jk)
        
        # For edge (j, k): add the cotangent from angle at vertex i.
        edge_key = tuple(sorted((j, k)))
        edge_weights[edge_key] = edge_weights.get(edge_key, 0.0) + cot_angle_at_i
        
        # For edge (i, k): add the cotangent from angle at vertex j.
        edge_key = tuple(sorted((i, k)))
        edge_weights[edge_key] = edge_weights.get(edge_key, 0.0) + cot_angle_at_j
        
        # For edge (i, j): add the cotangent from angle at vertex k.
        edge_key = tuple(sorted((i, j)))
        edge_weights[edge_key] = edge_weights.get(edge_key, 0.0) + cot_angle_at_k

    # Accumulate contributions for the Laplace operator at each vertex.
    # For each edge (i,j) with weight w, the contribution is w*(x_j - x_i) for vertex i,
    # and w*(x_i - x_j) for vertex j.
    for (i, j), w in edge_weights.items():
        laplace[i] += w * (points[j] - points[i])
        laplace[j] += w * (points[i] - points[j])
    
    # Compute the discrete Laplacian Δx and then mean curvature H = 0.5 * ||Δx||
    mean_curvature = np.zeros(N)
    for i in range(N):
        if vertex_area[i] > 1e-8:
            delta_x = laplace[i] / (2.0 * vertex_area[i])
            mean_curvature[i] = 0.5 * np.linalg.norm(delta_x)
        else:
            mean_curvature[i] = 0.0

    return mean_curvature

def process_single_file(args):
    """
    Process a single VTK file and save its mean curvature to H5
    """
    vtk_file, progress_queue = args
    try:
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Starting {vtk_file}")
        
        h5_file = os.path.splitext(vtk_file)[0] + '.h5'
        
        # Check if dataset already exists
        if os.path.exists(h5_file):
            with h5py.File(h5_file, 'r') as f:
                if 'mean_curvature' in f:
                    progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Skipped {vtk_file} (already processed)")
                    return f"Skipped {vtk_file}"

        # Process the file
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Reading mesh from {vtk_file}")
        mesh = read_vtk_file(vtk_file)
        
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Computing curvature for {vtk_file}")
        H = compute_mean_curvature(mesh)
        
        # Save to H5 file
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Saving results for {vtk_file}")
        with h5py.File(h5_file, 'a') as f:
            if 'mean_curvature' not in f:
                f.create_dataset('mean_curvature', data=H)
        
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] Completed {vtk_file}")
        return f"Successfully processed {vtk_file}"
        
    except Exception as e:
        error_msg = f"Error processing {vtk_file}: {str(e)}"
        progress_queue.put(f"[{datetime.now().strftime('%H:%M:%S')}] {error_msg}")
        return error_msg

def progress_monitor(queue):
    """Monitor and print progress messages from the queue"""
    while True:
        message = queue.get()
        if message == "DONE":
            break
        print(message, flush=True)

if __name__ == '__main__':
    # Add command line argument parsing
    parser = argparse.ArgumentParser(description='Compute mean curvature for VTK files')
    parser.add_argument('--serial', action='store_true', help='Run in serial mode instead of parallel')
    args = parser.parse_args()

    # Get all VTK files
    vtk_files = glob.glob('*.vtk')
    
    if not vtk_files:
        print("No VTK files found in current directory")
        sys.exit(0)
    
    if args.serial:
        print(f"Processing {len(vtk_files)} files in serial mode")
        # Process files serially
        for vtk_file in vtk_files:
            try:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] Processing {vtk_file}")
                h5_file = os.path.splitext(vtk_file)[0] + '.h5'
                
                # Check if dataset already exists
                if os.path.exists(h5_file):
                    with h5py.File(h5_file, 'r') as f:
                        if 'mean_curvature' in f:
                            print(f"[{datetime.now().strftime('%H:%M:%S')}] Skipped {vtk_file} (already processed)")
                            continue
                
                mesh = read_vtk_file(vtk_file)
                print(f"[{datetime.now().strftime('%H:%M:%S')}] Computing curvature")
                H = compute_mean_curvature(mesh)
                
                with h5py.File(h5_file, 'a') as f:
                    if 'mean_curvature' not in f:
                        f.create_dataset('mean_curvature', data=H)
                print(f"[{datetime.now().strftime('%H:%M:%S')}] Completed {vtk_file}")
                
            except Exception as e:
                print(f"Error processing {vtk_file}: {str(e)}")
    
    else:
        # Original parallel processing code
        num_processes = max(1, cpu_count() - 1)
        print(f"Processing {len(vtk_files)} files using {num_processes} processes")
        
        manager = Manager()
        progress_queue = manager.Queue()
        
        try:
            with Pool(processes=num_processes) as pool:
                # Start the progress monitor in the main process
                results = pool.map_async(process_single_file, 
                                       [(f, progress_queue) for f in vtk_files])
                
                # Monitor progress while waiting for results
                while not results.ready():
                    # Print any messages in the queue
                    while not progress_queue.empty():
                        print(progress_queue.get(), flush=True)
                    results.wait(timeout=1)
                
                # Get final results
                final_results = results.get()
                
                # Drain any remaining messages
                while not progress_queue.empty():
                    print(progress_queue.get(), flush=True)
                
        except KeyboardInterrupt:
            print("\nProcessing interrupted by user")
            sys.exit(1)
    
    print("\nAll files processed. Summary:")
    if not args.serial and 'final_results' in locals():
        for result in final_results:
            print(result)