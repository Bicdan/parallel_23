#pragma once

#include <vector>

class Handler {
    private:
        std::vector<std::vector<double>> _a, _b, _B;
        int M, N;
        double y_step, x_step, eps;

        double calculate_a(double y1, double y2, double x);
        double calculate_b(double y, double x1, double x2);
        double calculate_intersection_area(double x1, double x2, double y1, double y2);
        void init();

    public:
        Handler(int M, int N);
        const std::vector<std::vector<double>> getAw(const std::vector<std::vector<double>>& w);
        const std::vector<std::vector<double>> getB();
        void zeroing_borders(std::vector<std::vector<double>>& w);
};