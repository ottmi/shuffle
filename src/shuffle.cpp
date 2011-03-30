#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <sstream>
#include "globals.h"
#include "Alignment.h"

using namespace std;

int verbose = 0;

void parseArguments(int argc, char** argv, Options *options)
{
	char c;

	options->dataType = -1;
	options->randomizations = 100;
	options->minCo = 0;
	options->minPOC = 0.0;
	options->maxSmin = INT_MAX;
	options->maxEntropy = DBL_MAX;
	options->help = 0;
	options->grouping.push_back(0);

	while ( (c = getopt(argc, argv, "i:adg:r:s:o:c:p:m:e:v")) != -1)
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
			case 'g':
			{
				options->grouping.clear();
				int i;
				stringstream ss(optarg);
				while (ss >> i)
				{
					options->grouping.push_back(i);
					if (ss.peek() == ',')
						ss.ignore();
				}
				break;
			}
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
				options->minCo = atof(optarg);
				break;
			case 'p':
				options->minPOC = atof(optarg);
				break;
			case 'm':
				options->maxSmin = atoi(optarg);
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

	vector<int>::iterator min = min_element(options->grouping.begin(), options->grouping.end());
	vector<int>::iterator max = max_element(options->grouping.begin(), options->grouping.end());
	options->groupOffset = *min;
	options->groupLength = *max-*min+1;
}

void printSyntax()
{
	cout << "Syntax:" << endl;
	cout << "  shuffle -i <FILE> <-a|-d> [-g LIST] [-r NUM] [-o FILE [-c NUM] [-p NUM] [-m NUM] -[e NUM]] [-s FILE] [-v]" << endl;
	cout << "  shuffle -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -i\tInput alignment" << endl;
	cout << "  -a\tInput alignment is AA data" << endl;
	cout << "  -d\tInput alignment is DNA data" << endl;
	cout << "  -g\tGrouping of columns into sites, e.g. 0,1 for duplets and 0,1,2 for codons" << endl;
	cout << "  -r\tNumber of randomizations for POC computations [default: 100]" << endl;
	cout << "  -o\tOutput alignment" << endl;
	cout << "  -c\tMinimum Co score [default: 0]" << endl;
	cout << "  -p\tMinimum POC score [default: 0.0]" << endl;
	cout << "  -m\tMaximum Smin [default: " << INT_MAX << "]" << endl;
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

	Alignment alignment(&options);
	alignment.computeCompatibilityScores(options.randomizations);
	if (options.summaryFile.length())
		alignment.writeSummary(options.summaryFile);
	if (options.outputAlignment.length())
	{
		Alignment modifiedAlignment = alignment.getModifiedAlignment(options.minCo, options.minPOC, options.maxSmin, options.maxEntropy);
		modifiedAlignment.write(options.outputAlignment);
	}

	return 0;
}
