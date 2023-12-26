#pragma once

#include "Oracul.h"
#include "Info.h"
#include <mpi.h>

class MethodOfLeastRes {
public:
    static Info execute(Oracul& oracul, int worker_id, int num_workers, int M, int N, double delta = 1e-6);
};