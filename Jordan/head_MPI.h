#ifndef _JORD_
#define _JORD_
#include "mpi.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define EPS 1e-15

#define M 15

int solve(int , char *, int , int , MPI_Comm);

struct min_data{
	double a;
	int n;
};	

#endif
