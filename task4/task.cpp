#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "/opt/homebrew/opt/libomp/include/omp.h"

using namespace std;
using Matrix = vector<vector<double> >;

bool is_equal_double(double a, double b) {
    return abs(a - b) < 1e-8;
}

class Oracul {
    private:
        double X_left = 0., X_right = 3.;
        double Y_up = 4, Y_down = 0.;
        int M, N;
        double y_step, x_step, eps;

        Matrix a_matrix, b_matrix, B_matrix;

        double clip_x(double x) {
            return max(X_left, min(X_right, x));
        }

        double clip_y(double y) {
            return max(Y_down, min(Y_up, y));
        }

        double get_a_xy(double y1, double y2, double x) {
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

        double get_a(int i, int j) {
            double y1 = clip_y(Y_down + y_step * (i - 1. / 2));
            double y2 = clip_y(Y_down + y_step * (i + 1. / 2));
            double x = clip_x(X_left + x_step * (j - 1. / 2));
            return get_a_xy(y1, y2, x);
        }
        
        double get_b_xy(double y, double x1, double x2) {
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

        double get_b(int i, int j) {
            double x1 = clip_x(X_left + x_step * (j - 1. / 2));
            double x2 = clip_x(X_left + x_step * (j + 1. / 2));
            double y = clip_y(Y_down + y_step * (i - 1. / 2));
            return get_b_xy(y, x1, x2);
        }

        int is_inside_area(double x, double y) {
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

        double get_square(double x1, double x2, double y1, double y2) {
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

        void initialize_AB() {
            a_matrix = Matrix(M + 1, std::vector<double>(N + 1));
            b_matrix = Matrix(M + 1, std::vector<double>(N + 1));
            int i, j;
            #pragma omp parallel for private(i) private(j)
            for (i = 0; i < M + 1; ++i) {
                for (j = 0; j < N + 1; ++j) {
                    a_matrix[i][j] = get_a(i, j);
                    b_matrix[i][j] = get_b(i, j);
                }
            }
            B_matrix = Matrix(M + 1, std::vector<double>(N + 1));
            #pragma omp parallel for private(i) private(j)
            for (i = 0; i < M + 1; ++i) {
                for (j = 0; j < N + 1; ++j) {
                    double x1 = clip_x(X_left + x_step * (j - 1. / 2));
                    double x2 = clip_x(X_left + x_step * (j + 1. / 2));
                    double y1 = clip_y(Y_down + y_step * (i - 1. / 2));
                    double y2 = clip_y(Y_down + y_step * (i + 1. / 2));
                    B_matrix[i][j] = get_square(x1, x2, y1, y2) / x_step / y_step;
                }
            }
        }

    public:
        Oracul(int M, int N): M(M), N(N) {
            y_step = (Y_up - Y_down) / M;
            x_step = (X_right - X_left) / N;
            eps = max(x_step, y_step);
            eps *= eps;
            initialize_AB();
       }

       const Matrix getAw(Matrix &w) {
            Matrix Aw(M + 1, std::vector<double>(N + 1));
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

        const Matrix getB() {
            return B_matrix;
        }

        void border_conditions(Matrix &w) {
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
};

class Info {
    public:
        Matrix w, Aw, r;
        double r_norm_max;
    Info() {}
    Info(const Matrix w, const Matrix Aw, const Matrix r, double r_norm_max):
        w(w), Aw(Aw), r(r), r_norm_max(r_norm_max) {}
};

Info MethodOfLeastRes(Oracul &oracul, int M, int N, double delta = 1e-6) {
    int M_p = M + 1, N_p = N + 1;
    Matrix w(M_p, std::vector<double>(N_p));
    Info info;
    int max_iter = 1e6;
    double current_upd = 0.;

    for(int iter = 0; iter < max_iter; ++iter) {
        const Matrix& Aw = oracul.getAw(w);
        const Matrix& B = oracul.getB();
        Matrix r(M_p, std::vector<double>(N_p));
        int i, j;
        #pragma omp parallel for private(i) private(j)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                r[i][j] = Aw[i][j] - B[i][j];
            }
        }

        double tau = 0.;
        Matrix Ar = oracul.getAw(r);
        #pragma omp parallel for private(i) private(j) reduction(+: tau)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                tau += Ar[i][j] * r[i][j];
            }
        }

        double norm = 0.;
        #pragma omp parallel for private(i) private(j) reduction(+: norm)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                norm += Ar[i][j] * Ar[i][j];
            }
        }
        tau /= norm;

        #pragma omp parallel for private(i) private(j)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                w[i][j] -= tau * r[i][j];
            }
        }
        oracul.border_conditions(w);

        double r_norm_max = 0.;
        #pragma omp parallel for private(i) private(j) reduction(max: r_norm_max)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                r_norm_max = max(abs(r[i][j]), r_norm_max);
            }
        }

        current_upd = abs(tau) * r_norm_max;
            if (current_upd < delta || iter + 1 == max_iter) {
                info = Info(w, Aw, r, r_norm_max);
                break;
            }
    }
    //cout << current_upd << endl;
    return info;
}

void PrintMatrix(Matrix &w, string file_name = "output.txt") {
    std::ofstream fout(file_name);
    fout.precision(7);
    for (int i = w.size() - 1; i >= 0; --i) {
        for (int j = 0; j < w[0].size(); ++j) {
            fout << std::fixed << w[i][j] << " ";
        }
        fout << std::endl;
    }
    fout << std::endl;
    fout.close();
}

double run_test(int M, int N) {
    double min_time = 3600;
    int min_cycle = 5;
    for (int i = 1; i <= min_cycle; ++i) {
        double start = omp_get_wtime();
        Oracul oracul(M, N);
        Info info = MethodOfLeastRes(oracul, M, N);
        double end = omp_get_wtime();
        min_time = min(min_time, end - start);
        cout.precision(3);
        cout << "Time: " << std::fixed << end - start << endl;
        PrintMatrix(info.w);
    }
    return min_time;
}

int main() {
    int M = 20, N = 20;
    int num_threads = 20;
    // int sizes[] = {10, 20, 40};
    // int threads[] = {2, 4, 8, 16, 32};
    // num_threads = 32;
    // omp_set_num_threads(num_threads);
    // cout.precision(3);
    // for (int i = 0; i < 3; ++i) {
    //     M = sizes[i];
    //     N = sizes[i];
    //     double min_time = run_test(M, N);
    //     cout << "NumThreads: " << num_threads << "  M: " << M <<
    //         "  Time: " << std::fixed << min_time << endl;
    // }

    double min_time = run_test(M, N);
    cout << "NumThreads: " << num_threads << "  M: " << M <<
        "  Time: " << std::fixed << min_time << endl;

    return 0;
}