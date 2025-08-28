#include "helpers.h"
#include <random>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>

Point3D generateRandomPoint(std::mt19937& rng, std::uniform_real_distribution<double>& dist) {
    return {dist(rng), dist(rng), dist(rng)};
}

std::vector<Point3D> generateRandomPointsVector(int k, double min_coord, double max_coord) {
    // Seed the random number generator using the current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed); // Mersenne Twister engine

    // Define a uniform distribution for coordinates between min_coord and max_coord
    std::uniform_real_distribution<double> dist(min_coord, max_coord);

    std::vector<Point3D> random_points;
    random_points.reserve(k); // Pre-allocate memory for efficiency

    for (int i = 0; i < k; i++) {
        random_points.push_back(generateRandomPoint(rng, dist));
    }
    return random_points;
}

// Function to read points from a CSV-like file
std::vector<Point3D> read_points(const std::string& filename) {
    std::ifstream in(filename);
    std::vector<Point3D> points;
    std::string line;

    while (getline(in, line)) {
        if (line.empty()) continue; // skip blank lines

        std::stringstream ss(line);
        std::string token;
        std::vector<double> coords;

        // split by commas
        while (getline(ss, token, ',')) {
            coords.push_back(stod(token));
        }

        if (coords.size() == 3) {
            points.push_back({coords[0], coords[1], coords[2]});
        } else {
            std::cout << "Warning: line does not have 3 coords -> " << line << std::endl;
        }
    }
    return points;
}

void save_points(std::vector<Point3D> points, const std::string& filename) {
    std::ofstream points_stream(filename);
    if (points_stream.is_open()) {
        for (const auto& p : points) {
            points_stream << p.x << "," << p.y << "," << p.z << std::endl;
        }
    }
    else { std::cout << "Could not open file for writing" << std::endl; }
}

Point2D generateRandomUnitVector(std::mt19937& rng, std::normal_distribution<double>& dist) {
    double angle = dist(rng);
    return {cos(angle), sin(angle)};
}

std::vector<Point2D> generateRandomUnitVectors2D(int k, double offset, double std) {
    // Seed the random number generator using the current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed); // Mersenne Twister engine

    // Define a uniform distribution for coordinates between min_coord and max_coord
    // std::uniform_real_distribution<double> dist(0.0, 2*M_PI - EPS);
    std::normal_distribution<double> dist(offset, std);

    std::vector<Point2D> random_points;
    random_points.reserve(k); // Pre-allocate memory for efficiency

    for (int i = 0; i < k; ++i) {
        random_points.push_back(generateRandomUnitVector(rng, dist));
    }
    return random_points;
}
