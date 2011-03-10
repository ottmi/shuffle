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
	string outputAlignment;
	int dataType;
	int randomizations;
	string summaryFile;
	double minCompatibility;
	double minPOC;
	int maxMNIC;
	double maxEntropy;
} Options;

void parseArguments(int argc, char** argv, Options *options)
{
	char c;

	options->dataType = _DNA_DATA;
	options->randomizations = 100;
	options->minCompatibility = 0;
	options->minPOC = 0.0;
	options->maxMNIC = INT_MAX;
	options->maxEntropy = DBL_MAX;

	while ( (c = getopt(argc, argv, "adi:r:s:o:c:p:m:e:v")) != -1)
	{
		switch (c)
		{
			case 'a':
				options->dataType = _AA_DATA;
				break;
			case 'd':
				options->dataType = _DNA_DATA;
				break;
			case 'i':
				options->inputAlignment = optarg;
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
			default:
				cerr << "Unknown parameter: " << c << endl;
		}
	}
}

void printSyntax()
{
	exit(254);
}

int main(int argc, char** argv) {
	Options options;

	cout << PROGNAME << " " << VERSION << " [" << PROGDATE << "]" << endl;

	parseArguments(argc, argv, &options);
	if (!options.inputAlignment.length())
		printSyntax();

	Alignment alignment(options.inputAlignment, options.dataType);
	alignment.computeCompatibilityScores(options.randomizations);
	if (options.summaryFile.length())
		alignment.writeSummary(options.summaryFile);
	if (options.outputAlignment.length())
	{
		Alignment modifiedAlignment = alignment.getModifiedAlignment(options.minCompatibility, options.minPOC, options.maxMNIC, options.maxEntropy);
		modifiedAlignment.write("outfile.fasta");
	}

	return 0;
}
