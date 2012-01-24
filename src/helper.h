#ifndef HELPER_H_
#define HELPER_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
using namespace std;

int getMyId();
int getNumOfCpus();
void getMyShare(unsigned int &start, unsigned int &end, unsigned int n);
void createRandomSeed(unsigned short *seed, unsigned int n);

double factorial(int n);
string printTime(long t);
istream& safeGetline(istream& is, string& t);
string adjustString(string s, bool upercase=false);

#endif /* HELPER_H_ */
