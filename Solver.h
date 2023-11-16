#pragma once

#include "Handler.h"
#include "Printer.h"

class Solver {
public:
    static Printer execute(Handler& handler, int M, int N, double delta = 1e-6);
};