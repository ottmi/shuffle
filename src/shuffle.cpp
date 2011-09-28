#ifdef _OPENMP
#include <omp.h>
#endif

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

int parseArguments(int argc, char** argv, Options *options)
{
	char c;

	options->dataType = -1;
	options->removeDuplicates = 0;
	options->randomizations = 100;
	options->minCo = 0;
	options->minPOC = 0.0;
	options->maxSmin = INT_MAX;
	options->maxEntropy = DBL_MAX;
	options->hasMinMax = false;
	options->help = 0;
	options->grouping.push_back(0);
	options->windowSize = -1;
	options->windowStep = -1;

	int minGroup = 0;
	int maxGroup = 0;

	while ( (c = getopt(argc, argv, "i:t:g:d::b:w:r:s:o:c:p:m:e:n:v::h")) != -1)
	{
		switch (c)
		{
			case 'i':
				options->inputAlignment = optarg;
				break;
			case 't':
			{
				char type = optarg[0];
				switch (type)
				{
					case 'a':
					case 'A':
						options->dataType = _AA_DATA;
						break;
					case 'd':
					case 'D':
						options->dataType = _DNA_DATA;
						break;
					case 'n':
					case 'N':
						options->dataType = _ALPHANUM_DATA;
						break;
					default:
						cerr << "Unknown data type: " << optarg << endl;
						return 2;
						break;
				}
				break;
			}
			case 'g':
			{
				options->grouping.clear();
				minGroup = 255;
				int i;
				stringstream ss(optarg);
				while (ss >> i)
				{
					if (i >= 0)
						options->grouping.push_back(i);
					else
						i = -i;
					if (i > maxGroup)
						maxGroup = i;
					if (i < minGroup)
						minGroup = i;

					if (ss.peek() == ',')
						ss.ignore();
				}
				break;
			}
			case 'd':
				options->removeDuplicates = 1;
				if (optarg)
					options->reducedAlignment = optarg;
				break;
			case 'b':
				options->bowkersTest = optarg;
				break;
			case 'w':
			{
				int i;
				stringstream ss(optarg);
				while (ss >> i)
				{
					if (options->windowSize < 0)
						options->windowSize = i;
					options->windowStep = i;
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
				options->hasMinMax = true;
				break;
			case 'p':
				options->minPOC = atof(optarg);
				options->hasMinMax = true;
				break;
			case 'm':
				options->maxSmin = atoi(optarg);
				options->hasMinMax = true;
				break;
			case 'e':
				options->maxEntropy = atof(optarg);
				options->hasMinMax = true;
				break;
#ifdef _OPENMP
			case 'n':
				omp_set_num_threads(atoi(optarg));
				break;
#endif
			case 'v':
				if (optarg)
					verbose = atoi(optarg);
				else
					verbose = 1;
				break;
			case 'h':
				options->help = 1;
				break;
			default:
				if (c != '?')
					cerr << "Unknown parameter: " << c << endl;
				return 1;
		}
	}

	options->groupOffset = minGroup;
	options->groupLength = maxGroup - minGroup + 1;

	return 0;
}

void printSyntax()
{
	cout << "Syntax:" << endl;
	cout << "  shuffle -i<FILE> [-t<a|d|n>] [-g<LIST>] [-b<FILE> [-w<NUM>[,NUM]]] [-r<NUM>]" << endl;
	cout << "          [-o<FILE> [-c<NUM>] [-p<NUM>] [-m<NUM>] [-e<NUM>]]" << endl;
    cout << "          [-s<FILE>] [-d[FILE]] ";
#ifdef _OPENMP
	cout << "[-n<NUM>] ";
#endif
	cout << "[-v[NUM]]" << endl;

	cout << "  shuffle -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -i\tInput alignment" << endl;
	cout << "  -t\tData type a=AA, d=DNA, n=Alphanumeric [default: auto-detect]" << endl;
	cout << "  -g\tGrouping of sites, e.g. 0,1,-2 for duplets and 0,1,2 for codons" << endl;
	cout << "  -b\tPerform Bowker's test for symmetry and write results to file" << endl;
	cout << "  -w\tWindow size and step width for Bowker's test [default: sequence length]" << endl;
	cout << "  -r\tNumber of randomizations for POC computations [default: 100]" << endl;
	cout << "  -o\tOutput alignment" << endl;
	cout << "  -c\tMinimum Co score [default: 0]" << endl;
	cout << "  -p\tMinimum POC score [default: 0.0]" << endl;
	cout << "  -m\tMaximum Smin [default: " << INT_MAX << "]" << endl;
	cout << "  -e\tMaximum entropy [default: " << DBL_MAX << "]" << endl;
	cout << "  -s\tOutput file for site summary" << endl;
	cout << "  -d\tIdentify duplicates and optionally dump reduced alignment to file" << endl;
#ifdef _OPENMP
	cout << "  -n\tNumber of threads [default: " << omp_get_max_threads() << "]" << endl;
#endif
	cout << "  -v\tBe increasingly verbose" << endl;
	cout << "  -h\tThis help page" << endl;
	cout << endl;

	cout << "Note:" << endl;
	cout << "  .phy or .phylip extensions indicate sequential Phylip file format" << endl;
	cout << "  .fsa or .fasta extensions indicate Fasta file format" << endl;


	exit(254);
}

int main(int argc, char** argv) {
	Options options;

	cout << PROGNAME << " " << VERSION << "|";
#ifdef _OPENMP
	cout << "OpenMP|";
#endif
	cout << PROGDATE << endl << endl;;


	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;
	if (!options.inputAlignment.length() || options.help)
		printSyntax();

#ifdef _OPENMP
	cout << "Parallel execution with " << omp_get_max_threads() << " threads." << endl << endl;
#endif

	Alignment alignment(&options);

	if (options.removeDuplicates)
		alignment.removeDuplicates();

	if (options.reducedAlignment.length())
		alignment.write(options.reducedAlignment);

	if (options.summaryFile.length() || options.outputAlignment.length() || options.bowkersTest.length())
	{
		alignment.collectSites(&options);

		if (options.bowkersTest.length())
			alignment.computeBowkers(options.bowkersTest, options.windowSize, options.windowStep);

		if (options.summaryFile.length() || options.outputAlignment.length())
			alignment.collectInformativeSites(&options);

		if (options.summaryFile.length() || options.hasMinMax)
			alignment.computeCompatibilityScores(options.randomizations);

		if (options.summaryFile.length())
			alignment.writeSummary(options.summaryFile);

		if (options.outputAlignment.length())
		{
			Alignment modifiedAlignment = alignment.getModifiedAlignment(options.minCo, options.minPOC, options.maxSmin, options.maxEntropy);
			modifiedAlignment.write(options.outputAlignment);
		}
	}

	return 0;
}
