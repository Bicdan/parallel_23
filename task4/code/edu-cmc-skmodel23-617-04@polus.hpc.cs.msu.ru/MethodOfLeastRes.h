#pragma once

#include "Oracul.h"
#include "Info.h"

class MethodOfLeastRes {
public:
    static Info execute(Oracul& oracul, int M, int N, double delta = 1e-6);
};