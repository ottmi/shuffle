#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string>
#include <vector>
using namespace std;

#define PROGNAME "shuffle++"
#define VERSION "0.7.10"
#define PROGDATE "2011-12-01"

#define _DNA_DATA				0
#define	_AA_DATA				1
#define	_ALPHANUM_DATA			2

#define	_FASTA_FORMAT			0
#define	_PHYLIP_FORMAT			1

#define _DNA_MAP      "ACGTURYKMSWBDHVN?-"
#define _AA_MAP       "ACDEFGHIKLMNPQRSTVWYBJXZ?-"
#define _ALPHANUM_MAP "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?-"

extern int verbose;

typedef struct opt_struct
{
	string inputAlignment;
	int alignmentFormat;
	int dataType;
	string prefix;
	vector<int> grouping;
	int groupOffset;
	int groupLength;
	bool removeDuplicates;
	bool removeInformativeSitesDuplicates;
	bool convertAlignment;
	bool symmetryTest;
	bool writeExtendedTestResults;
	int windowSize;
	int windowStep;
	int randomizations;
	bool requireInformative;
	bool writeSiteSummary;
	bool writeRandomizedCo;
	bool filterAlignment;
	double minCo;
	double minPOC;
	int maxSmin;
	double maxEntropy;
	int help;
} Options;

#endif /* GLOBALS_H_ */
