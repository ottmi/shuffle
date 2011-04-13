#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <vector>
using namespace std;

#define PROGNAME "shuffle++"
#define VERSION  "0.30"
#define PROGDATE "2011-04-13"

extern int verbose;
extern int debug;

#define _DNA_DATA				0
#define	_AA_DATA				1
#define	_ALPHANUM_DATA	2

typedef struct opt_struct
{
	string inputAlignment;
	int dataType;
	string outputAlignment;
	int randomizations;
	string summaryFile;
	double minCo;
	double minPOC;
	int maxSmin;
	double maxEntropy;
	int help;
	vector<int> grouping;
	int groupOffset;
	int groupLength;
} Options;

#endif /* GLOBALS_H_ */
