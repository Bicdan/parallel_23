#include "Handler.h"
#include "Solver.h"
#include <iostream>
#include "/opt/homebrew/opt/libomp/include/omp.h"

void run_once() {
    int M = 40, N = 40;
    int num_threads = 20;
    omp_set_num_threads(num_threads);
    std::cout.precision(3);

    double time = 60.0 * 120;
    int num_cycles = 5;
    for (int i = 1; i <= num_cycles; ++i) {
        double start = omp_get_wtime();
        Handler handler(M, N);
        Solver solver = Solver();
        Printer info = solver.execute(handler, M, N);
        double end = omp_get_wtime();
        time = std::min(time, end - start);
        std::cout.precision(3);
        std::cout << "Time: " << std::fixed << end - start << std::endl;
        info.printMatrix(info.w, "output.txt");
    }
    std::cout << "M = " << M << ", N = " << N << ", NumThreads = " << num_threads << ". Time: " << std::fixed << time << std::endl;
}

double timeit(int M, int N) {
    double time = 60.0 * 120;
    int num_cycles = 5;
    for (int i = 1; i <= num_cycles; ++i) {
        double start = omp_get_wtime();
        Handler handler(M, N);
        Printer info = Solver().execute(handler, M, N);
        double end = omp_get_wtime();
        time = std::min(time, end - start);
    }
    return time;
}

void run_experiment() {
    int sizes[] = {10, 20, 40, 80};
    int threads[] = {2, 4, 8, 16, 32};
    int num_threads = 16;
    omp_set_num_threads(num_threads);
    std::cout.precision(3);
    for (int j = 0; j < 4; ++j) {
        int M = sizes[j];
        int N = sizes[j];
        double time = timeit(M, N);
        std::cout << "M = " << M << ", N = " << N << ", NumThreads = " << num_threads << ". Time: " << std::fixed << time << std::endl;
    }
}

int main() {
    //run_once();
    run_experiment();
    return 0;
}