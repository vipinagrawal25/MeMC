// Helper function to calculate polar angle for rotation
double polar_angle(const Vec3d& vector) {
    return atan2(sqrt(vector.x * vector.x + vector.y * vector.y), vector.z);
}

// Helper function to rotate a point around an axis by given angle
Vec3d rotate_point(const Vec3d& point, const Vec3d& axis, double theta) {
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    
    // Rodrigues' rotation formula
    Vec3d rotated = point * cos_theta + 
                   cross(axis, point) * sin_theta + 
                   axis * dot(axis, point) * (1.0 - cos_theta);
    return rotated;
}

// Helper function to sort 2D points by polar angle (anticlockwise)
vector<int> sort_2D_points_by_angle(const vector<Vec3d>& points) {
    vector<pair<double, int>> angle_index_pairs;
    
    for (int i = 0; i < points.size(); ++i) {
        double angle = atan2(points[i].y, points[i].x);
        angle_index_pairs.push_back({angle, i});
    }
    
    // Sort by angle
    sort(angle_index_pairs.begin(), angle_index_pairs.end());
    
    vector<int> sorted_indices;
    for (const auto& pair : angle_index_pairs) {
        sorted_indices.push_back(pair.second);
    }
    
    return sorted_indices;
}

/**
 * @brief Sort neighbors of each node in anticlockwise order by rotating coordinate system
 * @param R Array of position vectors for all nodes
 * @param Np Number of nodes
 * @param cmlst Cumulative list for neighbor indexing (cmlst[i] to cmlst[i+1] gives neighbors of node i)
 * @param node_nbr Array containing neighbor indices
 */
void sort_nbrs(Vec3d* R, int Np, int* cmlst, int* node_nbr) {
    Vec3d zhat(0.0, 0.0, 1.0);
    
    for (int i = 0; i < Np; ++i) {
        int start_idx = cmlst[i];
        int end_idx = cmlst[i + 1];
        int num_neighbors = end_idx - start_idx;
        
        if (num_neighbors <= 1) continue; // No need to sort if 0 or 1 neighbors
        
        Vec3d vector = R[i];
        
        // Calculate rotation axis (cross product of vector with z-axis)
        Vec3d vhat = cross(vector, zhat);
        double vnorm = magnitude(vhat);
        
        // If the vector is already lying at z-axis then there is no need to rotate
        if (vnorm > 1e-16) {
            vhat = vhat / vnorm;
            double theta = polar_angle(vector);
            
            // Rotate all the neighbors of this point
            vector<Vec3d> rotated_neighbors;
            vector<int> original_nbr_indices;
            
            for (int j = start_idx; j < end_idx; ++j) {
                int nbr_idx = node_nbr[j];
                Vec3d rotated = rotate_point(R[nbr_idx], vhat, theta);
                rotated_neighbors.push_back(rotated);
                original_nbr_indices.push_back(nbr_idx);
            }
            
            // Sort them in anticlockwise direction
            vector<int> sorted_indices = sort_2D_points_by_angle(rotated_neighbors);
            
            // Update node_nbr with sorted neighbors
            for (int j = 0; j < num_neighbors; ++j) {
                node_nbr[start_idx + j] = original_nbr_indices[sorted_indices[j]];
            }
        }
    }
}
