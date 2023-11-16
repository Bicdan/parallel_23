#pragma once

#include <vector>

class Handler {
    private:
        std::vector<std::vector<double>> _a, _b, _B;
        int M, N;
        double y_step, x_step, eps;

        double clip_x(double x);
        double clip_y(double y);
        double calculate_a(double y1, double y2, double x);
        double get_a(int i, int j);
        double calculate_b(double y, double x1, double x2);
        double get_b(int i, int j);
        int is_inside_area(double x, double y);
        double get_square(double x1, double x2, double y1, double y2);
        void init();

    public:
        Handler(int M, int N);
        const std::vector<std::vector<double>> getAw(const std::vector<std::vector<double>>& w);
        const std::vector<std::vector<double>> getB();
        void zeroing_borders(std::vector<std::vector<double>>& w);
};