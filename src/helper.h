#ifndef HELPER_H_
#define HELPER_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
using namespace std;

int getMyId();
int getNumOfCpus();

double factorial(int n);
string printTime(long t);
istream& safeGetline(istream& is, string& t);
string adjustString(string s, bool upercase=false);

#endif /* HELPER_H_ */
