#include <iostream>
#include <math.h>
using namespace std;
/*
From Numerical Recipes
*/

#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define FPMIN 1.0e-30
#define ITMAX 100

void gcf(double *gammcf, double a, double x, double *gln)
{
   double gammln(double xx);
   void nrerror(char error_text[]);
   int i;
   double an,b,c,d,del,h;

   *gln=gammln(a);
   b=x+1.0-a;
   c=1.0/FPMIN;
   d=1.0/b;
   h=d;
   for (i=1;i<=ITMAX;i++) {
      an = -i*(i-a);
      b += 2.0;
      d=an*d+b;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=b+an/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      if (fabs(del-1.0) < EPS) break;
   }
   if (i > ITMAX) cerr << "a too large, ITMAX too small in gcf" << endl;
   *gammcf=exp(-x+a*log(x)-(*gln))*h;
}


/*
From Numerical Recipes
*/

void gser(double *gamser, double a, double x, double *gln)
{
   double gammln(double xx);
   void nrerror(char error_text[]);
   int n;
   double sum,del,ap;

   *gln=gammln(a);
   if (x <= 0.0) {
      if (x < 0.0) cerr << "x less than 0 in routine gser" << endl;
      *gamser=0.0;
      return;
   } else {
      ap=a;
      del=sum=1.0/a;
      for (n=1;n<=ITMAX;n++) {
         ++ap;
         del *= x/ap;
         sum += del;
         if (fabs(del) < fabs(sum)*EPS) {
            *gamser=sum*exp(-x+a*log(x)-(*gln));
            return;
         }
      }
      cerr << "a too large, ITMAX too small in routine gser" << endl;
      return;
   }
}


/*
From Numerical Recipes
*/

double gammln(double xx)
{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
   int j;

   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.00000000019001;
   for (j=0;j<=5;j++) ser += cof[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}


/*
From Numerical Recipes
*/

double gammq(double a, double x)
{
   void gcf(double *gammcf, double a, double x, double *gln);
   void gser(double *gamser, double a, double x, double *gln);
   void nrerror(char error_text[]);
   double gamser,gammcf,gln;

   if (x < 0.0 || a <= 0.0) cerr << "Invalid arguments in routine gammq" << endl;
   if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return 1.0-gamser;
   } else {
      gcf(&gammcf,a,x,&gln);
      return gammcf;
   }
}
