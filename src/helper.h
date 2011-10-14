#ifndef HELPER_H_
#define HELPER_H_

#ifdef _GMP
#include <gmp.h>
#include <gmpxx.h>
#define bignum mpz_class
#else
#define bignum double
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include <string>
using namespace std;

bignum factorial(int n);
double bignum_log(bignum);
string printTime(long t);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);
double gammq(double a, double x);

#endif /* HELPER_H_ */
