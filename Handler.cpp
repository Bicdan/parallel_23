#include "Handler.h"
#include <algorithm>
#include <iostream>

Handler::Handler(int M, int N) : M(M), N(N) {
    y_step = (Y_up - Y_down) / M;
    x_step = (X_right - X_left) / N;
    eps = std::max(x_step, y_step);
    eps *= eps;
    initialize_AB();
}

bool is_equal_double(double a, double b, double epsilon = 1e-8) {
    return std::abs(a - b) < epsilon;
}

double Handler::clip_x(double x) {
    return std::max(X_left, std::min(X_right, x));
}

double Handler::clip_y(double y) {
    return std::max(Y_down, std::min(Y_up, y));
}

double Handler::get_a_xy(double y1, double y2, double x) {
    if (y1 > -4*x/3 + 4) {
        return 1 / eps;
    }
    if (y2 <= -4*x/3 + 4) {
        return 1;
    } else {
        double y = -4*x / 3 + 4;
        double l = y - y1;
        return l / y_step+ (1 - l / y_step) / eps;
    }
}

double Handler::get_a(int i, int j) {
    double y1 = clip_y(Y_down + y_step * (i - 1. / 2));
    double y2 = clip_y(Y_down + y_step * (i + 1. / 2));
    double x = clip_x(X_left + x_step * (j - 1. / 2));
    return get_a_xy(y1, y2, x);
}

double Handler::get_b_xy(double y, double x1, double x2) {
    if (x1 > -3*y/4 + 3) {
        return 1 / eps;
    }
    if (x2 <= -3*y/4 + 3) {
        return 1;
    } else {
        double x = -3*y / 4 + 3;
        double l = x - x1;
        return l / y_step+ (1 - l / y_step) / eps;
    }
}

double Handler::get_b(int i, int j) {
    double x1 = clip_x(X_left + x_step * (j - 1. / 2));
    double x2 = clip_x(X_left + x_step * (j + 1. / 2));
    double y = clip_y(Y_down + y_step * (i - 1. / 2));
    return get_b_xy(y, x1, x2);
}

int Handler::is_inside_area(double x, double y) {
    if (x > 0 && y > 0 && 3*y + 4*x < 12) {
        return 1;
    } else {
        if (is_equal_double(x, 0) && y >= 0 && y <= 4) {
            return 0;
        } else if (is_equal_double(y, 0) && x >= 0 && x <= 3) {
            return 0;
        } else if (x > 0 && y > 0 && is_equal_double(3*y + 4*x, 12)) {
            return 0;
        }
    }
    return -1;
}

double Handler::get_square(double x1, double x2, double y1, double y2) {
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

void Handler::initialize_AB() {
    a_matrix = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));
    b_matrix = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));

    #pragma omp parallel for
    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            a_matrix[i][j] = get_a(i, j);
            b_matrix[i][j] = get_b(i, j);
        }
    }

    B_matrix = std::vector<std::vector<double>>(M + 1, std::vector<double>(N + 1));

    #pragma omp parallel for
    for (int i = 0; i < M + 1; ++i) {
        for (int j = 0; j < N + 1; ++j) {
            double x1 = clip_x(X_left + x_step * (j - 1. / 2));
            double x2 = clip_x(X_left + x_step * (j + 1. / 2));
            double y1 = clip_y(Y_down + y_step * (i - 1. / 2));
            double y2 = clip_y(Y_down + y_step * (i + 1. / 2));
            B_matrix[i][j] = get_square(x1, x2, y1, y2) / x_step / y_step;
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
                double a_ij = a_matrix[i][j], b_ij = b_matrix[i][j];
                double a_i1j = a_matrix[i + 1][j];
                double b_ij1 = b_matrix[i][j + 1];
                Aw[i][j] = -1. / y_step *
                    (a_i1j / y_step * (w[i + 1][j] - w[i][j]) -
                        a_ij / y_step * (w[i][j] - w[i - 1][j])) +
                    - 1. / x_step *
                    (b_ij1 / x_step * (w[i][j + 1] - w[i][j]) -
                        b_ij / x_step * (w[i][j] - w[i][j - 1]));
            }
        }
    }
    return Aw;
}

const std::vector<std::vector<double>> Handler::getB() {
    return B_matrix;
}

void Handler::border_conditions(std::vector<std::vector<double>>& w) {
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
