#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include "globals.h"
#include "Alignment.h"

using namespace std;

int verbose = 0;

typedef struct opt_struct
{
	string inputAlignment;
	int dataType;
	string outputAlignment;
	int randomizations;
	string summaryFile;
	double minCompatibility;
	double minPOC;
	int maxMNIC;
	double maxEntropy;
	int help;
} Options;

void parseArguments(int argc, char** argv, Options *options)
{
	char c;

	options->dataType = -1;
	options->randomizations = 100;
	options->minCompatibility = 0;
	options->minPOC = 0.0;
	options->maxMNIC = INT_MAX;
	options->maxEntropy = DBL_MAX;
	options->help = 0;

	while ( (c = getopt(argc, argv, "adi:r:s:o:c:p:m:e:v")) != -1)
	{
		switch (c)
		{
			case 'i':
				options->inputAlignment = optarg;
				break;
			case 'a':
				options->dataType = _AA_DATA;
				break;
			case 'd':
				options->dataType = _DNA_DATA;
				break;
			case 'r':
				options->randomizations = atoi(optarg);
				break;
			case 's':
				options->summaryFile = optarg;
				break;
			case 'o':
				options->outputAlignment = optarg;
				break;
			case 'c':
				options->minCompatibility = atof(optarg);
				break;
			case 'p':
				options->minPOC = atof(optarg);
				break;
			case 'm':
				options->maxMNIC = atoi(optarg);
				break;
			case 'e':
				options->maxEntropy = atof(optarg);
				break;
			case 'v':
				verbose = 1;
				break;
			case 'h':
				options->help = 1;
				break;
			default:
				cerr << "Unknown parameter: " << c << endl;
		}
	}
}

void printSyntax()
{
	cout << "Syntax:" << endl;
	cout << "  shuffle -i <FILE> <-a|-d> [-r NUM] [-o FILE [-c NUM] [-p NUM] [-m NUM] -[e NUM]] [-s FILE] [-v]" << endl;
	cout << "  shuffle -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -i\tInput alignment" << endl;
	cout << "  -a\tInput alignment is AA data" << endl;
	cout << "  -d\tInput alignment is DNA data" << endl;
	cout << "  -r\tNumber of randomizations for POC computations [default: 100]" << endl;
	cout << "  -o\tOutput alignment" << endl;
	cout << "  -c\tMinimum compatibility score [default: 0]" << endl;
	cout << "  -p\tMinimum POC score [default: 0.0]" << endl;
	cout << "  -m\tMaximum MNIC [default: " << INT_MAX << "]" << endl;
	cout << "  -e\tMaximum entropy [default: " << DBL_MAX << "]" << endl;
	cout << "  -s\tOutput file for site summary" << endl;
	cout << "  -v\tBe verbose " << endl;
	cout << "  -h\tThis help page" << endl;
	cout << endl;

	cout << "Note:" << endl;
	cout << "  .phy or .phylip extensions indicate sequential Phylip file format" << endl;
	cout << "  .fsa or .fasta extensions indicate Fasta file format" << endl;


	exit(254);
}

int main(int argc, char** argv) {
	Options options;

	cout << PROGNAME << " " << VERSION << " [" << PROGDATE << "]" << endl << endl;

	parseArguments(argc, argv, &options);
	if (!(options.inputAlignment.length() && options.dataType >= 0) || options.help)
		printSyntax();

	Alignment alignment(options.inputAlignment, options.dataType);
	alignment.computeCompatibilityScores(options.randomizations);
	if (options.summaryFile.length())
		alignment.writeSummary(options.summaryFile);
	if (options.outputAlignment.length())
	{
		Alignment modifiedAlignment = alignment.getModifiedAlignment(options.minCompatibility, options.minPOC, options.maxMNIC, options.maxEntropy);
		modifiedAlignment.write(options.outputAlignment);
	}

	return 0;
}
