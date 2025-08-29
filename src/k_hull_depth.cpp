#include <Rcpp.h>
using namespace Rcpp;
#include "helpers.h"
#include "load_points.h"
#include <cassert>

long long SDk(const std::vector<Point3D>& Xinput, const Point3D& x, int k) {
    int n = Xinput.size();
    std::vector<Point3D> X;
    X.reserve(n);
    for (const auto& p : Xinput) {
        Point3D q{p.x - x.x, p.y - x.y, p.z - x.z};
        if (fabs(q.x) > EPS) X.push_back(q);
    }
    n = (int)X.size();
    if (n == 0 || k > n) return 0LL;

    struct HelperStruct {
        double x, y;
        int color, label;
        double sort_key;
    };

    struct PointY {
        double y1, y2;
        int sign;
    };
    std::vector<PointY> Y; Y.reserve(n);
    int n_up = 0;
    for (const auto& p : X) {
        double y1 = p.x / p.z;
        double y2 = p.y / p.z;
        int s = sgn(p.z);
        Y.push_back({y1, y2, s});
        if (s == 1) n_up++;
    }
    int n_down = n - n_up;

    long long SDk = combinations(n, k) - combinations(n_up, k) - combinations(n_down, k);
    if (n == 1) return SDk;

    std::vector<long long> counts(n - 1, 0LL);

    std::vector<HelperStruct> sorted_points(n-1);

    for (int idx = 0; idx < n; idx++) {
        int j = 0;
        for (int i = 0; i < n; i++) {
            if (i != idx) {
                HelperStruct* sorted_point = &sorted_points[j];
                PointY* pointy_i = &Y[i], *pointy_idx = &Y[idx];
                sorted_point->x = pointy_i->y1 - pointy_idx->y1,
                sorted_point->y = pointy_i->y2 - pointy_idx->y2,
                sorted_point->color = (pointy_i->sign != pointy_idx->sign) ? 1 : 0,
                sorted_point->label = (sorted_point->y * pointy_i->sign > 0.0) ? 1 : 0,
                sorted_point->sort_key = -(sorted_point->x / sorted_point->y);
                j++;
            }
        }
        int m = (int)sorted_points.size();
        if (m == 0) continue;

        std::sort(sorted_points.begin(), sorted_points.end(), [](const auto& a, const auto& b) {
            return a.sort_key < b.sort_key;
        });

        int EL = 0;
        for (int j = 1; j < m; j++) EL += sorted_points[j].label;
        counts[EL] += sorted_points[0].color;

        for (int j = 0; j < m - 1; j++) {
            if (sorted_points[j].label == 0) EL++;
            if (sorted_points[j+1].label == 1) EL--;
            counts[EL] += sorted_points[j+1].color;
        }
    }

    long long SUM = 0LL;
    for (int l = 1; l <= n - 1; l++) {
        long long cnt = counts[l - 1];
        if (cnt == 0) continue;
        SUM += cnt * (combinations(l - 1, k - 2) + combinations(n - l - 1, k - 2));
    }

    SDk -= SUM / 4;
    return SDk;
}

// [[Rcpp::export]]
double SDk_wrapper(NumericMatrix Xinput_mat, NumericVector x_vec, int k) {
    int n = Xinput_mat.nrow();
    if (Xinput_mat.ncol() != 3 || x_vec.size() != 3) {
        stop("Xinput must be n x 3 and x must be length 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> Xinput(n);
    for (int i = 0; i < n; i++) {
        Xinput[i].x = Xinput_mat(i, 0);
        Xinput[i].y = Xinput_mat(i, 1);
        Xinput[i].z = Xinput_mat(i, 2);
    }

    // Convert vector to Point3D
    Point3D x = { x_vec[0], x_vec[1], x_vec[2] };

    long long result = SDk(Xinput, x, k);

    return static_cast<double>(result) / combinations(n, k); // R does not have long long
}

long long SDk_parallel(const std::vector<Point3D>& Xinput, const Point3D& x, int k) {
    int n = Xinput.size();
    std::vector<Point3D> X;
    X.reserve(n);
    for (const auto& p : Xinput) {
        Point3D q{p.x - x.x, p.y - x.y, p.z - x.z};
        if (fabs(q.x) > EPS) X.push_back(q);
    }
    n = (int)X.size();
    if (n == 0 || k > n) return 0LL;

    struct HelperStruct {
        double x, y;
        int color, label;
        double sort_key;
    };

    struct PointY {
        double y1, y2;
        int sign;
    };
    std::vector<PointY> Y; Y.reserve(n);
    int n_up = 0;
    for (const auto& p : X) {
        double y1 = p.x / p.z;
        double y2 = p.y / p.z;
        int s = sgn(p.z);
        Y.push_back({y1, y2, s});
        if (s == 1) n_up++;
    }
    int n_down = n - n_up;

    long long SDk = combinations(n, k) - combinations(n_up, k) - combinations(n_down, k);
    if (n == 1) return SDk;

    std::vector<long long> counts(n - 1, 0LL);
    std::vector<HelperStruct> sorted_points(n-1);
#pragma omp parallel
    {
    std::vector<long long> local_counts(n-1, 0LL);  // private to each thread
    std::vector<HelperStruct> local_sorted_points(n-1);

    #pragma omp for schedule(dynamic)
    for (int idx = 0; idx < n; idx++) {
        int j = 0;
        for (int i = 0; i < n; i++) {
            if (i != idx) {
                HelperStruct* sorted_point = &local_sorted_points[j];
                PointY* pointy_i = &Y[i], *pointy_idx = &Y[idx];
                sorted_point->x = pointy_i->y1 - pointy_idx->y1;
                sorted_point->y = pointy_i->y2 - pointy_idx->y2;
                sorted_point->color = (pointy_i->sign != pointy_idx->sign) ? 1 : 0;
                sorted_point->label = (sorted_point->y * pointy_i->sign > 0.0) ? 1 : 0;
                sorted_point->sort_key = -(sorted_point->x / sorted_point->y);
                j++;
            }
        }
        // int m = (int)local_sorted_points.size();
        int m = n-1;
        if (m == 0) continue;

        std::sort(local_sorted_points.begin(), local_sorted_points.end(),
                  [](const auto& a, const auto& b) {
                      return a.sort_key < b.sort_key;
                  });

        int EL = 0;
        for (int j = 1; j < m; j++) EL += local_sorted_points[j].label;
        local_counts[EL] += local_sorted_points[0].color;

        for (int j = 0; j < m - 1; j++) {
            if (local_sorted_points[j].label == 0) EL++;
            if (local_sorted_points[j+1].label == 1) EL--;
            local_counts[EL] += local_sorted_points[j+1].color;
        }
    }

    // Merge local_counts into global counts
    #pragma omp critical
    {
        for (int i = 0; i < n - 1; i++)
            counts[i] += local_counts[i];
    }
    }
    long long SUM = 0LL;
    for (int l = 1; l <= n - 1; l++) {
        long long cnt = counts[l - 1];
        if (cnt == 0) continue;
        SUM += cnt * (combinations(l - 1, k - 2) + combinations(n - l - 1, k - 2));
    }

    SDk -= SUM / 4;
    return SDk;
}

// [[Rcpp::export]]
double SDk_parallel_wrapper(NumericMatrix Xinput_mat, NumericVector q, int k) {
    int n = Xinput_mat.nrow();
    if (Xinput_mat.ncol() != 3 || q.size() != 3) {
        stop("Xinput must be n x 3 and x must be length 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> Xinput(n);
    for (int i = 0; i < n; i++) {
        Xinput[i].x = Xinput_mat(i, 0);
        Xinput[i].y = Xinput_mat(i, 1);
        Xinput[i].z = Xinput_mat(i, 2);
    }

    // Convert vector to Point3D
    Point3D x = { q[0], q[1], q[2] };

    long long result = SDk_parallel(Xinput, x, k);

    return static_cast<double>(result) / combinations(n, k); // R does not have long long
}
