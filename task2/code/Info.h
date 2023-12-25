#pragma once

#include <vector>
#include <fstream>

class Info {
public:
    std::vector<std::vector<double>> w, Aw, r;
    double r_norm_max;

    Info();
    Info(const std::vector<std::vector<double>>& w, const std::vector<std::vector<double>>& Aw,
         const std::vector<std::vector<double>>& r, double r_norm_max);

    static void printMatrix(const std::vector<std::vector<double>>& w, const std::string& file_name) {
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
};