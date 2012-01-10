#ifndef HELPER_H_
#define HELPER_H_

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include <string>
using namespace std;

double factorial(int n);
string printTime(long t);
istream& safeGetline(istream& is, string& t);
string adjustString(string s, bool upercase=false);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);
double gammq(double a, double x);

#endif /* HELPER_H_ */
