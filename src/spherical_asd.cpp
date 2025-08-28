#include <Rcpp.h>
#include "helpers.h"
#include "load_points.h"
#include "simplicial_depth_2d.h"
#include "simplicial_depth_bruteforce.h"
#include <cmath>
#include <cassert>
using namespace Rcpp;

// The main algorithm for fast 3D angular simplicial depth
long long spherical_asd(Point3D ray, const std::vector<Point3D>& P) {
    int n = P.size();
    if (n < 3) return 0;

    // Step 1: Project points to a unit sphere.
    std::vector<Point3D> P_prime;
    P_prime.reserve(n);
    for (const auto& pi : P) {
        if (magnitude_sq(pi) > EPS) {
            P_prime.push_back(normalize(pi));
        }
    }
    n = P_prime.size();
    if (n < 3) return 0;

    bool ray_color = (ray.z > 0);
    // Project ray onto the z=1 plane
    Point2D ray_projected = {ray.x / ray.z, ray.y / ray.z};
    // Project all the other points onto this plane
    std::vector<Point2D> P_projected(n);
    for (int i=0; i<n; i++) P_projected[i] = {P[i].x / P[i].z, P[i].y / P[i].z};
    // Find the color of each point (is it in the same half space as ray?)
    std::vector<bool> P_colors(n);
    for (int i=0; i<n; i++) P_colors[i] = ((P[i].z > 0) == ray_color);
    // Shift all the points so that the projected ray becomes the origin
    for (int i=0; i<n; i++) P_projected[i] = P_projected[i] - ray_projected;

    std::vector<Point2D> P_red;
    std::vector<Point2D> P_blue;
    P_red.reserve(n);
    P_blue.reserve(n);
    for (int i=0; i<n; i++) {
        if (P_colors[i] == 0) P_blue.push_back(P_projected[i]);
        else P_red.push_back(P_projected[i]);
    }
    int n_blue_points = P_blue.size();
    int n_red_points = P_red.size();

    
    long long N1 = simplicial_depth_2d({0.0, 0.0}, P_red); // Number of triangles of the opposite color of ray that contain the origin

    // Helper struct for the rotational sweep algorithm
    struct ProjPoint {
        Point2D p_2d;
        double slope; // Angle on the projection plane
        bool label;
        int index;
    };

    long long N2 = 0;
    #pragma omp for schedule(dynamic)
    for (const auto& red_point : P_red) {
        Point2D center_shifted = {-red_point.x, -red_point.y};
        std::vector<ProjPoint> P_blue_shifted(n_blue_points);
        bool reference_label = center_shifted.y > 0 ? 1 : 0;
        double reference_slope = -center_shifted.x / center_shifted.y;
        for (int j=0; j<n_blue_points; j++) {
            P_blue_shifted[j].p_2d = P_blue[j] - red_point;
            P_blue_shifted[j].label = ((P_blue_shifted[j].p_2d.y > 0) == reference_label);
            P_blue_shifted[j].slope = -P_blue_shifted[j].p_2d.x / P_blue_shifted[j].p_2d.y;
            P_blue_shifted[j].index = j;
        }

        std::sort(P_blue_shifted.begin(), P_blue_shifted.end(), [](const ProjPoint& a, const ProjPoint& b){return a.slope < b.slope;});

        int reference_order = 0;
        while (reference_order < n_blue_points) {
            if (P_blue_shifted[reference_order].slope < reference_slope) reference_order += 1;
            else break;
        }

        // Case 1: 0 *1* 0
        // Count zeros to the left of reference order, to the right of ref. order and then multiply them
        int s_left = 0;
        int s_right = 0;
        for (int j=0; j<reference_order; j++) {
            if (P_blue_shifted[j].label == 0) s_left++;
        }
        for (int j=reference_order; j<n_blue_points; j++) {
            if (P_blue_shifted[j].label == 0) s_right++;
        }
        long long case_1 = s_left * s_right;


        // Case 2: 1 0 *1*
        // Counting number of 1 0 subsequences to the left of reference_order
        long long case_2 = 0, ones = 0;
        for (int j=0; j<reference_order; j++) {
            if (P_blue_shifted[j].label == 1) ones++;
            else case_2 += ones;
        }

        // Case 2: *1* 0 1
        // Counting number of 1 0 subsequences to the left of reference_order
        long long case_3 = 0, zeros = 0;
        for (int j=reference_order; j<n_blue_points; j++) {
            if (P_blue_shifted[j].label == 0) zeros++;
            else case_3 += zeros;
        }

        #pragma omp critical
        {
        N2 = N2 + case_1 + case_2 + case_3;
        }
    }

    // Number of intersections between all red segments with blue segments whose one endpoint is the origin
    long long N3 = 0;
    #pragma omp for schedule(dynamic)
    for (const auto& blue_point : P_blue) {
        // Rotate all red points counterclockwise so that the blue point Blue[b,] ends up on the negative x-axis
        double theta = atan2(blue_point.y, blue_point.x);
        double alpha = M_PI - theta;
        Matrix2D rotation_matrix = {{cos(alpha), -sin(alpha)}, {sin(alpha), cos(alpha)}};

        std::vector<ProjPoint> P_red_rotated(n_red_points);
        int red_labels_ones = 0, red_labels_zeros = 0;
        for (int j=0; j<n_red_points; j++) {
            P_red_rotated[j].p_2d = rotation_matrix * P_red[j];
            P_red_rotated[j].slope = -P_red_rotated[j].p_2d.x / P_red_rotated[j].p_2d.y;
            P_red_rotated[j].label = (P_red_rotated[j].p_2d.y > 0);
            red_labels_ones += (P_red_rotated[j].label == 1);
            red_labels_zeros += (P_red_rotated[j].label == 0);
        }
        Point2D blue_point_rotated = rotation_matrix * blue_point;
        std::sort(P_red_rotated.begin(), P_red_rotated.end(), [](const ProjPoint& a, const ProjPoint& b) {return a.slope < b.slope;});
        
        long long total_right = 0, ones = 0;
        for (const auto& red_point : P_red_rotated) {
            if (red_point.label == 1) ones++;
            else total_right += ones;
        }

        for (int j=0; j<n_red_points; j++) {
            P_red_rotated[j].p_2d = P_red_rotated[j].p_2d - blue_point_rotated;
            P_red_rotated[j].slope = -P_red_rotated[j].p_2d.x / P_red_rotated[j].p_2d.y;
        }
        std::sort(P_red_rotated.begin(), P_red_rotated.end(), [](const ProjPoint& a, const ProjPoint& b) {return a.slope < b.slope;});

        long long total_left = 0, zeros = 0;
        for (const auto& red_point : P_red_rotated) {
            if (red_point.label == 0) zeros++;
            else total_left += zeros;
        }
        #pragma omp critical
        {
        N3 = N3 + red_labels_ones * red_labels_zeros - total_right - total_left;
        }
    }
    return N1 + N2 + N3;
}

// [[Rcpp::export]]
double spherical_asd_wrapper(NumericVector ray_vec, NumericMatrix P_mat) {
    int n = P_mat.nrow();
    if (P_mat.ncol() != 3 || ray_vec.size() != 3) {
        stop("P must be n x 3 and ray must be length 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> P(n);
    for (int i = 0; i < n; i++) {
        P[i].x = P_mat(i, 0);
        P[i].y = P_mat(i, 1);
        P[i].z = P_mat(i, 2);
    }

    // Convert vector to Point3D
    Point3D ray = { ray_vec[0], ray_vec[1], ray_vec[2] };

    long long result = spherical_asd(ray, P);

    return static_cast<double>(result);
}

// int main() {
//
//     Point3D ray = {0.0, 0.0, 1.0};
//     // Point3D ray = {-0.5267046,  0.7868826, -0.3215558};
//
//     std::vector<Point3D> points3d = generateRandomPointsVector(500, (double)-10.0, (double)10.0);
//     // for (int i=0; i<points3d.size(); i++) points3d[i] = normalize(points3d[i]);
//     // std::vector<Point3D> points3d = {{1.0, 1.0, 1.0},
//     //                                  {1.0, -1.0, 1.0},
//     //                                  {-1.0, -0.25, 1.0},
//     //                                  {-0.5, 0.25, -1.0}};
//     // std::vector<Point3D> points3d = {{6.0, 2.85243, 0.00378952}, {-6.64376, 4.91298, -6.13063}, {6.92123, 1.92036, -8.65979}, {2.52833, 7.83858, -3.2764}, {-6.36378, -8.72196, 7.22423}, {8.35234, 0.265809, -0.0881092}};
//     // std::vector<Point3D> points3d = {{-3.33763, 5.0594, 9.92002}, {2.54541, 4.50546, 3.8714}, {-7.6481, -9.06328, -7.00999}, {8.06797, -6.40542, -1.14897}};
//     // std::vector<Point3D> points3d = {{-7.29813, -1.22975, 4.72286}, {5.48261, -7.92464, -0.667692}, {9.47525, 3.51214, 9.50944}, {9.27595, 5.16611, -8.48985}};
//     // std::vector<Point3D> points3d = {{0.479484, -0.793659, 0.374434}, {-0.493119, -0.690162, 0.529632}, {-0.741655, 0.653489, 0.151331}, {-0.573914, -0.0269235, -0.818473}, {0.597051, -0.631851, -0.494261}, {0.325311, 0.137545, 0.93555}, {0.866613, -0.369835, -0.334969}, {-0.630406, 0.068966, 0.773196}, {0.283099, -0.643798, -0.7109}, {0.618698, -0.246906, 0.745822}};
//     // std::vector<Point3D> points3d = {{1.46989, -9.00436, -3.25317}, {9.1502, 0.448838, -2.74683}, {-7.76597, -3.8567, 8.63104}, {4.54706, -4.70863, -9.51804}, {6.6199, 6.90817, 4.26871}, {-7.05371, 0.00665808, 3.95459}, {-1.13463, 8.99006, -5.27225}, {-6.04753, 9.13125, 9.89727}, {-7.76234, 6.7066, -5.77162}, {-9.95967, 4.52668, 0.864219}};
//
//     auto slow_start = std::chrono::high_resolution_clock::now();
//     long long slow_triangles = spherical_asd_bruteforce(ray, points3d);
//     // long long slow_triangles = angular_simplicial_depth_3d_fast(ray, points3d);
//     // long long slow_triangles = 0;
//     std::cout << "SD3D - brute force: " << slow_triangles << std::endl;
//     auto slow_stop = std::chrono::high_resolution_clock::now();
//     auto fast_start = std::chrono::high_resolution_clock::now();
//     long long fast_triangles = spherical_asd(ray, points3d);
//     std::cout << "SD3D - fast algo:   " << fast_triangles << std::endl;
//     auto fast_stop = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> slow_time = slow_stop - slow_start;
//     std::chrono::duration<double, std::milli> fast_time = fast_stop - fast_start;
//     std::cout << "Slow time: " << slow_time.count() << " ms" << std::endl;
//     std::cout << "Fast time:  " << fast_time.count() << " ms" << std::endl;
//
//     // if (slow_triangles != fast_triangles) {
//     //     std::cout << "Points:" << std::endl << "{";
//     //     for (const auto& p : points3d) std::cout << p << ", ";
//     //     std::cout << "}" << std::endl;
//     // }
//
//     return 0;
// }
