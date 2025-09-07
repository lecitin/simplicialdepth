#include <Rcpp.h>
#include "helpers.h"
#include <cassert>
#include <random>
#include <chrono>
using namespace Rcpp;

// Helper struct for sorting 2D points by angle
struct Point2D_indexed {
    Point2D p;
    double angle;
    double angle_mod_PI;
    bool label;
};

// An O(n log(n)) algorithm for angular simplicial depth
// Very similar to the regular algorithm, but the sweeping line only goes until it meets the -ray
long long circular_asd(Point2D ray, const std::vector<Point2D>& points) {
    int m = points.size();
    if (m < 2) return 0;

    double ray_angle = atan2(ray.y, ray.x);
    ray_angle = ray_angle < 0 ? ray_angle + 2 * M_PI : ray_angle;

    std::vector<Point2D_indexed> sorted_2d_points(m);
    #pragma omp for schedule(static)
    for(int i=0; i<m; ++i) {
        sorted_2d_points[i].p = points[i] - ray;
        double angle = atan2(points[i].y, points[i].x);
        sorted_2d_points[i].angle = angle < 0 ? angle + 2 * M_PI : angle;
        sorted_2d_points[i].angle -= ray_angle;
        sorted_2d_points[i].angle = sorted_2d_points[i].angle < 0 ? sorted_2d_points[i].angle + 2 * M_PI : sorted_2d_points[i].angle;
    }

    double ray_antipodal_angle = atan2(-ray.y, -ray.x);
    ray_antipodal_angle = ray_antipodal_angle < 0 ? ray_antipodal_angle + 2 * M_PI : ray_antipodal_angle;

    // std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
    //     return a.angle < b.angle;
    // });
    parallel_sort(sorted_2d_points, [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle < b.angle;
    });

    long long non_containing_triangles = 0;
    int k = 1, i = 0; // Two-pointer for finding the antipodal point
    for (i = 0; i < m && sorted_2d_points[i].angle < M_PI; ++i) {
        if (k < i + 1) k = i + 1;
        while (k < i + m) {
            double angle_diff = sorted_2d_points[k % m].angle - sorted_2d_points[i].angle;
            if (angle_diff < 0) angle_diff += 2 * M_PI;
            if (angle_diff >= M_PI - EPS) break;
            k++;
        }
        int num_in_half_plane = k - i - 1;
        non_containing_triangles += combinations(num_in_half_plane, 1);
    }
    while (k < i + m + 10) {
        double angle_diff = sorted_2d_points[k % m].angle - M_PI;
        if (angle_diff < 0) angle_diff += 2 * M_PI;
        if (angle_diff >= M_PI - EPS) break;
        k++;
    }
    int num_in_half_plane = k - i;
    non_containing_triangles += combinations(num_in_half_plane, 2);

    return combinations(m, 2) - non_containing_triangles;
}

// An O(n log(n)) algorithm for angular simplicial depth based on the 01 sequence method
long long circular_asd_01(Point2D ray, const std::vector<Point2D>& points) {
    int m = points.size();
    if (m < 2) return 0;

    double ray_angle = atan2(ray.y, ray.x);
    ray_angle = ray_angle < 0 ? ray_angle + 2 * M_PI : ray_angle;

    int N0 = 0, N1 = 0;

    std::vector<Point2D_indexed> sorted_2d_points(m);
    #pragma omp for schedule(static)
    for(int i=0; i<m; ++i) {
        sorted_2d_points[i].p = points[i];
        // Find the polar angle with respect to the ray
        sorted_2d_points[i].angle = atan2pi(points[i]) - ray_angle;
        sorted_2d_points[i].angle = sorted_2d_points[i].angle < 0 ? sorted_2d_points[i].angle + 2 * M_PI : sorted_2d_points[i].angle;
        sorted_2d_points[i].label = (sorted_2d_points[i].angle < M_PI);
        if (sorted_2d_points[i].label == 0) N0++;
        else N1++;
        // Angles modulo PI
        sorted_2d_points[i].angle_mod_PI = sorted_2d_points[i].angle < M_PI ? sorted_2d_points[i].angle : sorted_2d_points[i].angle - M_PI;
    }

    // std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
    //     return a.angle_mod_PI < b.angle_mod_PI;
    // });
    parallel_sort(sorted_2d_points, [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle_mod_PI < b.angle_mod_PI;
    });

    long long subsequences_10 = 0;
    // Count the number of "10" subsequences
    for (const auto& p : sorted_2d_points) {
        if (p.label == 1) subsequences_10 += N0;
        else N0--;
    }
    return subsequences_10;
}

// [[Rcpp::export]]
long long circular_asd(NumericMatrix X, NumericVector ray) {
    int n = X.nrow();
    std::vector<Point2D> P(n);
    for (int i=0; i<n; i++) {
        P.push_back({X(i,0), X(i,1)});
    }
    Point2D xx = {ray[0], ray[1]};
    return circular_asd_01(xx, P);
}

// An O(n log(n)) algorithm for angular simplicial depth based on the 01 sequence method
std::vector<Arc> circular_asd_01_all_points(const std::vector<Point2D>& points) {
    int m = points.size();
    std::vector<Arc> arcs(m);
    if (m < 2) return arcs;
    std::vector<long long> N10s(m);
    std::vector<double> angles(m);

    int N0 = 0, N1 = 0;

    std::vector<Point2D_indexed> sorted_2d_points(m);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<m; i++) {
        sorted_2d_points[i].p = points[i],
        // Find the polar angle with respect to the ray
        sorted_2d_points[i].angle = atan2pi(points[i]),
        sorted_2d_points[i].label = (sorted_2d_points[i].angle < M_PI);
        // Angles modulo PI
        sorted_2d_points[i].angle_mod_PI = sorted_2d_points[i].angle < M_PI ? sorted_2d_points[i].angle : sorted_2d_points[i].angle - M_PI;
        angles[i] = sorted_2d_points[i].angle;
    }

    for (const auto& p : sorted_2d_points) N1 += p.label;
    N0 = m - N1;
    int N0_total = N0, N1_total = N1;

    // std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
    //     return a.angle_mod_PI < b.angle_mod_PI;
    // });
    parallel_sort(sorted_2d_points, [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle_mod_PI < b.angle_mod_PI;
    });

    // Compute the first N10
    long long N10 = 0;
    for (const auto& p : sorted_2d_points) {
        if (p.label == 1) N10 += N0;
        else N0--;
    }
    // Save as first entry of all points
    N10s[0] = N10;

    for (int i=0, j=1; i<2*m-1 && j<m; i++) {
        if (sorted_2d_points[mod(i,m)].label == i / m) {
            N0_total = N0_total - 1, N1_total = N1_total + 1;
        }
        else {
            N1_total = N1_total - 1;
            N10 = N10 - N0_total + N1_total;
            N0_total = N0_total + 1;
            N10s[j] = N10;
            j++;
        }
    }

    // Sort the angles by size (this only requires O(2*n)=O(n), since the presorted points only need to be iterated through twice
    for (int i=0, j=0; i<2*m && j<m; i++) {
        if (sorted_2d_points[mod(i,m)].label == 1-i/m)
            angles[j++] = sorted_2d_points[mod(i,m)].angle;
    }

    #pragma omp parallel for schedule(static)
    for (int i=0; i<m; i++) {
        arcs[i].left_angle = angles[mod(i-1,m)],
        arcs[i].right_angle = angles[i],
        arcs[i].depth = N10s[i];
    }
    
    return arcs;
}

// [[Rcpp::export]]
DataFrame circular_asd_all_arcs(NumericMatrix points_mat) {
    int n = points_mat.nrow();
    std::vector<Point2D> points_vec(n);

    // Convert matrix to std::vector<Point2D>
    for (int i = 0; i < n; i++) {
        points_vec[i].x = points_mat(i, 0);
        points_vec[i].y = points_mat(i, 1);
    }

    std::vector<Arc> arcs = circular_asd_01_all_points(points_vec);

    long long combs = combinations(n, 2);

    // Convert std::vector<Arc> to DataFrame
    int m = arcs.size();
    NumericVector start_angles(m), end_angles(m), depths(m);
    for (int i = 0; i < m; i++) {
        start_angles[i] = arcs[i].left_angle;
        end_angles[i]   = arcs[i].right_angle;
        depths[i]       = static_cast<double>(arcs[i].depth) / combs;
    }

    return DataFrame::create(
        Named("start") = start_angles,
        Named("end")   = end_angles,
        Named("depth") = depths
    );
}

//
// int main() {
//
//     unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//     std::mt19937 rng(seed); // Mersenne Twister engine
//
//     // Define a uniform distribution for coordinates between min_coord and max_coord
//     std::uniform_real_distribution<double> dist(0.0, 2*M_PI - EPS);
//
//     int m = 1000000;
//
//     // Generate a random ray (point where we want to calculate the ASD)
//     // double q_angle = dist(rng);
//     double q_angle = 0.0;
//     Point2D q = {cos(q_angle), sin(q_angle)};
//
//     // Generate the circular data
//     std::vector<Point2D> points2d = generateRandomUnitVectors2D(m, 4.0, 2.0);
//
//     // Time the brute force algorithm
//     auto brute_start = std::chrono::high_resolution_clock::now();
//     long long brute=0;
//     brute = circular_asd_01(q, points2d);
//     // for (int i=0; i<m-1; i++) circular_asd(q, points2d);
//     auto brute_stop = std::chrono::high_resolution_clock::now();
//
//     // Time the fast algorithm
//     auto fast_start = std::chrono::high_resolution_clock::now();
//     // std::vector<Arc> ASD_fast = circular_asd_01_all_points(points2d);
//     // long long fast = ASD[0].depth;
//     long long fast = circular_asd(q, points2d);
//     auto fast_stop = std::chrono::high_resolution_clock::now();
//
//     // Compare the two results
//     if (brute != fast) std::cout << "Brute: " << brute << std::endl << "Fast: " << fast << std::endl;
//     assert(brute == fast);
//     
//     // Average the computation time of each algorithm
//     std::chrono::duration<double, std::milli> brute_time = brute_stop - brute_start;
//     std::chrono::duration<double, std::milli> fast_time = fast_stop - fast_start;
//
//     std::cout << "Bruteforce time: " << brute_time.count() << "ms\n";
//     std::cout << "Fast algorithm time: " << fast_time.count() << "ms\n";
//
//     std::cout << std::endl;
//     // Find arcs with min and max depth
//     // Arc min_arc = *std::min_element(ASD_fast.begin(), ASD_fast.end(),
//     //     [](const Arc& a, const Arc& b){ return a.depth < b.depth; });
//     // Arc max_arc = *std::max_element(ASD_fast.begin(), ASD_fast.end(),
//     //     [](const Arc& a, const Arc& b){ return a.depth < b.depth; });
//
//
//     // Print the min and max depth arc midpoints
//     // std::cout << "Angle min: " << "(" << cos(min_arc.mid_point()) << ", " << sin(min_arc.mid_point()) << ")" << ", depth: " << double(min_arc.depth) / double(combinations(m,2)) << std::endl;
//     // std::cout << "Angle max: " << "(" << cos(max_arc.mid_point()) << ", " << sin(max_arc.mid_point()) << ")" << ", depth: " << double(max_arc.depth) / double(combinations(m,2)) << std::endl;
//     
//     return 0;
// }
