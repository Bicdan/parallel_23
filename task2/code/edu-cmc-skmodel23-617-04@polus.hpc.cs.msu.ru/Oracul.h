#pragma once

#include <vector>

class Oracul {
    private:
        double X_left = 0., X_right = 3.;
        double Y_up = 4, Y_down = 0.;
        int M, N;
        double y_step, x_step, eps;

        std::vector<std::vector<double>> a_matrix, b_matrix, B_matrix;

        double clip_x(double x);
        double clip_y(double y);
        double get_a_xy(double y1, double y2, double x);
        double get_a(int i, int j);
        double get_b_xy(double y, double x1, double x2);
        double get_b(int i, int j);
        int is_inside_area(double x, double y);
        double get_square(double x1, double x2, double y1, double y2);
        void initialize_AB();

    public:
        Oracul(int M, int N);
        const std::vector<std::vector<double>> getAw(const std::vector<std::vector<double>>& w);
        const std::vector<std::vector<double>> getB();
        void border_conditions(std::vector<std::vector<double>>& w);
};