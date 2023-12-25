#include "Info.h"
#include <fstream>

Info::Info() {}

Info::Info(const std::vector<std::vector<double>>& w, const std::vector<std::vector<double>>& Aw,
           const std::vector<std::vector<double>>& r, double r_norm_max)
    : w(w), Aw(Aw), r(r), r_norm_max(r_norm_max) {}