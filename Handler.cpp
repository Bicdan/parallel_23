#include "Handler.h"
#include <algorithm>
#include <iostream>

Handler::Handler(int M, int N) : M(M), N(N) {
    h1 = 3. / N;
    h2 = 4. / M;
    eps = std::max(h1, h2) * std::max(h1, h2);
    init();
}

double Handler::calculate_a(double y1, double y2, double x) {
    if (y1 > -4*x/3 + 4) {
        return 1 / eps;
    }
    if (y2 <= -4*x/3 + 4) {
        return 1;
    } else {
        double y = -4*x / 3 + 4;
        double l = y - y1;
        return l / h2 + (1 - l / h2) / eps;
    }
}

double Handler::calculate_b(double y, double x1, double x2) {
    if (x1 > -3*y/4 + 3) {
        return 1 / eps;
    }
    if (x2 <= -3*y/4 + 3) {
        return 1;
    } else {
        double x = -3*y / 4 + 3;
        double l = x - x1;
        return l / h1 + (1 - l / h1) / eps;
    }
}

double Handler::calculate_intersection_area(double x1, double x2, double y1, double y2) {
    if (3*y1 + 4*x1 > 12) {
        return 0;
    } else if (3*y2 + 4*x2 < 12) {
        return (x2 - x1) * (y2 - y1);
    } else {
        double x = -3*y1 / 4 + 4;
        double y = -4*x1 / 3 + 4;
        if (y <= y2 && x <= x2) {
            return (y - y1) * (x - x1) / 2;
        } else if (y > y2 && x <= x2) {
            double xx = -3*y2 / 4 + 4;
            return (y - y1) * (x - x1) / 2 - (y - y2) * (xx - x1) / 2;
        } else if (y <= y2 && x > x2) {
            double yy = -4*x2 / 3 + 4;
            return (y - y1) * (x - x1) / 2 - (yy - y1) * (x - x2) / 2;
        } else if (y > y2 && x > x2) {
            double xx = -3*y2 / 4 + 4;
            double yy = -4*x2 / 3 + 4;
            return (y - y1) * (x - x1) / 2 - (y - y2) * (xx - x1) / 2 - (yy - y1) * (x - x2) / 2;
        }
    }
    return (x2 - x1) * (y2 - y1) / 2;
}

void Handler::init() {
    _a = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));
    _b = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));

    #pragma omp parallel for
    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            double y1 = h2 * (i - 0.5);
            double y2 = h2 * (i + 0.5);
            double x = h1 * (j - 0.5);
            _a[i][j] = calculate_a(y1, y2, x);

            double x1 = h1 * (j - 0.5);
            double x2 = h1 * (j + 0.5);
            double y = h2 * (i - 0.5);
            _b[i][j] = calculate_b(y, x1, x2);
        }
    }

    _B = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));

    #pragma omp parallel for
    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            double x1 = h1 * (j - 0.5);
            double x2 = h1 * (j + 0.5);
            double y1 = h2 * (i - 0.5);
            double y2 = h2 * (i + 0.5);
            _B[i][j] = calculate_intersection_area(x1, x2, y1, y2) / h1 / h2;
        }
    }
}

const std::vector<std::vector<double>> Handler::getAw(const std::vector<std::vector<double>>& w) {
    std::vector<std::vector<double>> Aw(M + 1, std::vector<double>(N + 1));
    int i, j;
    #pragma omp parallel for private(i) private(j)
    for (i = 0; i < M + 1; ++i) {
        for (j = 0; j < N + 1; ++j) {
            if (i > 0 && j > 0 && i < M && j < N) {
                double a_ij = _a[i][j];
                double b_ij = _b[i][j];
                double a_i1j = _a[i + 1][j];
                double b_ij1 = _b[i][j + 1];
                Aw[i][j] = -1. / h2 *
                    (a_i1j / h2 * (w[i + 1][j] - w[i][j]) -
                        a_ij / h2 * (w[i][j] - w[i - 1][j])) +
                    - 1. / h1 *
                    (b_ij1 / h1 * (w[i][j + 1] - w[i][j]) -
                        b_ij / h1 * (w[i][j] - w[i][j - 1]));
            }
        }
    }
    return Aw;
}

const std::vector<std::vector<double>> Handler::getB() {
    return _B;
}

void Handler::zeroing_borders(std::vector<std::vector<double>>& w) {
    int i, j;
    #pragma omp parallel for private(i)
    for (i = 0; i < M + 1; ++i) {
        w[i][0] = 0;
        w[i][N] = 0;
    }
    #pragma omp parallel for private(j)
    for (j = 0; j < N + 1; ++j) {
        w[0][j] = 0;
        w[M][j] = 0;
    }
}
