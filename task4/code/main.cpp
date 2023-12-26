#include "Oracul.h"
#include "MethodOfLeastRes.h"
#include <iostream>
#include "/opt/homebrew/opt/libomp/include/omp.h"
#include <mpi.h>

int num_workers, worker_id;

int main() {
    int M = 20, N = 20;
    int num_threads = 20;
    std::cout.precision(3);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

    double min_time = 3600;
    int min_cycle = 5;
    omp_set_num_threads(num_threads);
    for (int i = 1; i <= min_cycle; ++i) {
        double start = omp_get_wtime();
        Oracul oracul(M, N);
        MethodOfLeastRes method = MethodOfLeastRes();
        Info info = method.execute(oracul, M, N);
        double end = omp_get_wtime();
        min_time = std::min(min_time, end - start);
        std::cout.precision(3);
        std::cout << "Time: " << std::fixed << end - start << std::endl;
        info.printMatrix(info.w, "output.txt");
    }
    std::cout << "NumThreads: " << num_threads << "  M: " << M << " N: " << N <<
        "  Time: " << std::fixed << min_time << std::endl;

    MPI_Finalize();
    return 0;
}