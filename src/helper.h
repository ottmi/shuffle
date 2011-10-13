/*
 * helper.h
 *
 *  Created on: 15/09/2011
 *      Author: ott029
 */

#ifndef HELPER_H_
#define HELPER_H_

#ifdef _GMP
#include <gmp.h>
#include <gmpxx.h>
#define bignum mpz_class
#else
#define bignum double
#endif

bignum factorial(int n);
double bignum_log(bignum);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammln(double xx);
double gammq(double a, double x);

#endif /* HELPER_H_ */
