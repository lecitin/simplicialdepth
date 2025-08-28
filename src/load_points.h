#include "helpers.h"
#include <random>

Point3D generateRandomPoint(std::mt19937& rng, std::uniform_real_distribution<double>& dist);
std::vector<Point3D> generateRandomPointsVector(int k, double min_coord, double max_coord);
std::vector<Point3D> read_points(const std::string& filename);
void save_points(const std::vector<Point3D> points, const std::string& filename);
Point2D generateRandomUnitVector(std::mt19937& rng, std::normal_distribution<double>& dist);
std::vector<Point2D> generateRandomUnitVectors2D(int k, double offset, double std);
