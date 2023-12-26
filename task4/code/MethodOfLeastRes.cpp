#include "MethodOfLeastRes.h"
#include <iostream>
#include "Oracul.h"

Info MethodOfLeastRes::execute(Oracul& oracul, int worker_id, int num_workers, int M, int N, double delta) {
    int grid_size = 1;
    while (grid_size * grid_size < num_workers) {
        ++grid_size;
    }

    MPI_Status status;
    
    int M_p = M + 1, N_p = N + 1;
    std::vector<std::vector<double>> w(M_p, std::vector<double>(N_p));
    Info info;
    int max_iter = 1e6;
    double current_upd = 0.;

    for(int iter = 0; iter < max_iter; ++iter) {
        broadcast(w, buf, &status, M, N, worker_id, grid_size, num_workers);
        
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
        MPI_Allreduce(&tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        std::vector<std::vector<double>> Ar = oracul.getAw(r);
        #pragma omp parallel for private(i) private(j) reduction(+: tau)
        for (i = 0; i < M_p; ++i) {
            for (j = 0; j < N_p; ++j) {
                tau += Ar[i][j] * r[i][j];
            }
        }

        double norm = 0.;
        MPI_Allreduce(&norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
        broadcast(r, buf, &status, M, N, worker_id, grid_size, num_workers);

        current_upd = abs(tau) * r_norm_max;
        if (current_upd < delta || iter + 1 == max_iter) {
            info = Info(w, Aw, r, r_norm_max);
            break;
        }
        MPI_Allreduce(&current_upd, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    return info;
}
