#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <vector>
using namespace std;

#define PROGNAME "shuffle++"
#define VERSION "0.6.6"
#define PROGDATE "2011-09-30"
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
	string prefix;
	vector<int> grouping;
	int groupOffset;
	int groupLength;
	bool removeDuplicates;
	bool removeInformativeSitesDuplicates;
	bool symmetryTest;
	bool writeExtendedTestResults;
	int windowSize;
	int windowStep;
	int randomizations;
	bool writeSiteSummary;
	bool filterAlignment;
	double minCo;
	double minPOC;
	int maxSmin;
	double maxEntropy;
	int help;
} Options;

#endif /* GLOBALS_H_ */
