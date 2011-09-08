#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <vector>
using namespace std;

#define PROGNAME "shuffle++"
#define VERSION "0.5.2"
#define PROGDATE "2011-09-08"
extern int verbose;

#define _DNA_DATA				0
#define	_AA_DATA				1
#define	_ALPHANUM_DATA	2

#define _DNA_MAP      "ACGTURYKMSWBDHVN?-"
#define _AA_MAP       "ACDEFGHIKLMNPQRSTVWYBJXZ?-"
#define _ALPHANUM_MAP "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?-"

typedef struct opt_struct
{
	string inputAlignment;
	int dataType;
	int removeDuplicates;
	string reducedAlignment;
	string outputAlignment;
	int randomizations;
	string summaryFile;
	double minCo;
	double minPOC;
	int maxSmin;
	double maxEntropy;
	bool hasMinMax;
	int help;
	vector<int> grouping;
	int groupOffset;
	int groupLength;
} Options;

#endif /* GLOBALS_H_ */
