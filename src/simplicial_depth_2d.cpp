#include "helpers.h"

// An O(n log(n)) algorithm for regular simplicial depth
// [[Rcpp::export]]
long long simplicial_depth_2d(Point2D q, const std::vector<Point2D>& points) {
    int m = points.size();
    if (m < 3) return 0;

    struct Point2D_indexed {
        Point2D p;
        double angle;
    };

    std::vector<Point2D_indexed> sorted_2d_points(m);
    for(int i=0; i<m; i++) {
        sorted_2d_points[i].p = points[i] - q;
        sorted_2d_points[i].angle = atan2pi(points[i].y - q.y, points[i].x - q.x);
    }

    std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle < b.angle;
    });

    long long non_containing_triangles = 0;
    int k = 1; // Two-pointer for finding the antipodal point
    for (int i = 0; i < m; i++) {
        if (k < i + 1) k = i + 1;
        while (k < i + m) {
            double angle_diff = sorted_2d_points[k % m].angle - sorted_2d_points[i].angle;
            if (angle_diff < 0) angle_diff += 2 * M_PI;
            if (angle_diff >= M_PI - EPS) break;
            k++;
        }
        int num_in_half_plane = k - i - 1;
        non_containing_triangles += combinations(num_in_half_plane, 2);
    }

    return combinations(m, 3) - non_containing_triangles;
}
