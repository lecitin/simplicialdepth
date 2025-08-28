#include "helpers.h"

double signed_volume(Point3D p1, Point3D p2, Point3D p3, Point3D p4) {
    return (1.0 / 6.0) * dot(cross(p2 - p1, p3 - p1), p4 - p1);
}

bool is_inside_tetrahedron(Point3D p1, Point3D p2, Point3D p3, Point3D p4, Point3D q) {
    double c1 = signed_volume(p1, p2, p3, q);
    double c2 = signed_volume(p1, p2, q,  p4);
    double c3 = signed_volume(p1, q,  p3, p4);
    double c4 = signed_volume(q,  p2, p3, p4);

    bool has_neg = (c1 < 0) || (c2 < 0) || (c3 < 0) || (c4 < 0);
    bool has_pos = (c1 > 0) || (c2 > 0) || (c3 > 0) || (c4 > 0);
    
    // Inside if all cross products have the same sign (or zero)
    return !(has_neg && has_pos);
}

long long simplicial_depth_bruteforce(const Point3D& p, const std::vector<Point3D>& P) {
    int n = P.size();
    if (n < 4) return 0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                for (int l = k + 1; l < n; l++) {
                    if (is_inside_tetrahedron(P[i], P[j], P[k], P[l], p)) {
                        count++;
                    }
                }
            }
        }
    }
    return count;
}

double signed_area(Point2D p1, Point2D p2, Point2D p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

// Returns a boolean value indicating whether q lies inside a triangle formed by p1, p2 and p3
bool is_inside_triangle(Point2D p1, Point2D p2, Point2D p3, Point2D q) {
    double c1 = signed_area(p1, p2, q);
    double c2 = signed_area(p1, q,  p3);
    double c3 = signed_area(q , p2, p3);

    bool has_neg = (c1 < 0) || (c2 < 0) || (c3 < 0);
    bool has_pos = (c1 > 0) || (c2 > 0) || (c3 > 0);
    
    // Inside if all cross products have the same sign (or zero)
    return !(has_neg && has_pos);
}

// A straightforward O(n^3) implementation of retular SD for verification
long long simplicial_depth_bruteforce(const Point2D& p, const std::vector<Point2D>& P) {
    int n = P.size();
    if (n < 3) return 0;
    long long count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                if (is_inside_triangle(P[i], P[j], P[k], p)) {
                    count++;
                }
            }
        }
    }
    return count;
}

// A straightforward O(n^2) implementation of angular SD for verification of the fast ASD algorithm
long long circular_asd_bruteforce(const Point2D& ray, const std::vector<Point2D>& P) {
    Point2D p = {-ray.x, -ray.y};
    Point2D origin = {0.0, 0.0};
    int n = P.size();
    if (n < 2) return 0;
    long long count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (is_inside_triangle(P[i], P[j], p, origin)) {
                count++;
            }
        }
    }
    return count;
}

int spherical_asd_bruteforce(const Point3D& ray, const std::vector<Point3D>& P) {
    int n = P.size();
    Point3D ray_antipodal = {-ray.x, -ray.y, -ray.z};
    if (n < 4) return 0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            for (int k = j + 1; k < n; k++) {
                if (is_inside_tetrahedron(P[i], P[j], P[k], ray_antipodal, {0.0, 0.0, 0.0})) {
                    count++;
                }
            }
        }
    }
    return count;
}
