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
	options->help = 0;
	options->grouping.push_back(0);

	int minGroup = 0;
	int maxGroup = 0;

	while ( (c = getopt(argc, argv, "i:t:g:d::r:s:o:c:p:m:e:v::h")) != -1)
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
						options->dataType = _AA_DATA;
						break;
					case 'd':
						options->dataType = _DNA_DATA;
						break;
					case 'n':
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
	cout << "  shuffle -i<FILE> -t<a|d|n> [-d[FILE]] [-g<LIST>] [-r<NUM>] [-o<FILE> [-c<NUM>] [-p<NUM>] [-m<NUM>] -[e<NUM>]] [-s<FILE>] [-v[NUM]]" << endl;
	cout << "  shuffle -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -i\tInput alignment" << endl;
	cout << "  -t\tInput alignment data type a=AA, d=DNA, n=Alphanumeric" << endl;
	cout << "  -d\tRemove duplicates and optionally dump reduced alignment to file" << endl;
	cout << "  -g\tGrouping of columns into sites, e.g. 0,1 for duplets and 0,1,2 for codons" << endl;
	cout << "  -r\tNumber of randomizations for POC computations [default: 100]" << endl;
	cout << "  -o\tOutput alignment" << endl;
	cout << "  -c\tMinimum Co score [default: 0]" << endl;
	cout << "  -p\tMinimum POC score [default: 0.0]" << endl;
	cout << "  -m\tMaximum Smin [default: " << INT_MAX << "]" << endl;
	cout << "  -e\tMaximum entropy [default: " << DBL_MAX << "]" << endl;
	cout << "  -s\tOutput file for site summary" << endl;
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

	cout << PROGNAME << " " << VERSION << " [" << PROGDATE << "]" << endl << endl;

#ifdef _OPENMP
	cout << "This is the OpenMP version, running in parallel on " << omp_get_max_threads() << " Threads." << endl;
	cout << endl;
#endif

	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;
	if (!(options.inputAlignment.length() && options.dataType >= 0) || options.help)
		printSyntax();


	Alignment alignment(&options);

	if (options.removeDuplicates)
		alignment.removeDuplicates();

	if (options.reducedAlignment.length())
		alignment.write(options.reducedAlignment);

	alignment.collectInformativeSites(&options);

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
