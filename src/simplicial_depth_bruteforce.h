#include "helpers.h"

double signed_volume(Point3D p1, Point3D p2, Point3D p3, Point3D p4);
bool is_inside_tetrahedron(Point3D p1, Point3D p2, Point3D p3, Point3D p4, Point3D q);
long long simplicial_depth_bruteforce(const Point3D& p, const std::vector<Point3D>& P);
double signed_area(Point2D p1, Point2D p2, Point2D p3);
bool is_inside_triangle(Point2D p1, Point2D p2, Point2D p3, Point2D q);
long long circular_asd_bruteforce(const Point2D& ray, const std::vector<Point2D>& P);
long long simplicial_depth_bruteforce(const Point2D& p, const std::vector<Point2D>& P);
int spherical_asd_bruteforce(const Point3D& ray, const std::vector<Point3D>& P);
