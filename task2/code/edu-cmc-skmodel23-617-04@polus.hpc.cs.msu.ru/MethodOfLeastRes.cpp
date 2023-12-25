#include "MethodOfLeastRes.h"
#include <iostream>
#include "Oracul.h"

Info MethodOfLeastRes::execute(Oracul& oracul, int M, int N, double delta) {
    int M_p = M + 1, N_p = N + 1;
    std::vector<std::vector<double>> w(M_p, std::vector<double>(N_p));
    Info info;
    int max_iter = 1e6;
    double current_upd = 0.;

    for(int iter = 0; iter < max_iter; ++iter) {
        const std::vector<std::vector<double>>& Aw = oracul.getAw(w);
        const std::vector<std::vector<double>>& B = oracul.getB();
        std::vector<std::vector<double>> r(M_p, std::vector<double>(N_p));
        int i, j;
        #pragma omp parallel for private(i) private(j)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                r[i][j] = Aw[i][j] - B[i][j];
            }
        }

        double tau = 0.;
        std::vector<std::vector<double>> Ar = oracul.getAw(r);
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
                r_norm_max = std::max(abs(r[i][j]), r_norm_max);
            }
        }

        current_upd = abs(tau) * r_norm_max;
            if (current_upd < delta || iter + 1 == max_iter) {
                info = Info(w, Aw, r, r_norm_max);
                break;
            }
    }
    return info;
}
