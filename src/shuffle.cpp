#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <sstream>
#include <cstdio>
#include "globals.h"
#include "Alignment.h"
#include "helper.h"

using namespace std;

int verbose = 0;
int numProcs = 1;
int myId = 0;


int parseArguments(int argc, char** argv, Options *options)
{
	char c;

	if (argc < 2)
		return 0;

	options->inputAlignment = string(argv[--argc]);

	options->dataType = -1;
	options->removeDuplicates = false;
	options->removeInformativeSitesDuplicates = false;
	options->symmetryTest = false;
	options->writeExtendedTestResults = false;
	options->randomizations = 100;
	options->writeSiteSummary = false;
	options->filterAlignment = false;
	options->minCo = 0;
	options->minPOC = 0.0;
	options->maxSmin = INT_MAX;
	options->maxEntropy = DBL_MAX;
	options->help = 0;
	options->grouping.push_back(0);
	options->windowSize = -1;
	options->windowStep = -1;
	options->writeExtendedTestResults = false;

	int minGroup = 0;
	int maxGroup = 0;

	while ( (c = getopt(argc, argv, "t:p:g:diyxw:r:sf:n:v::h")) != -1)
	{
		switch (c)
		{
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
			case 'p':
				options->prefix = optarg;
				break;
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
				options->removeDuplicates = true;
				break;
			case 'i':
				options->removeInformativeSitesDuplicates = true;
				break;
			case 'y':
				options->symmetryTest = true;
				break;
			case 'x':
				options->writeExtendedTestResults = true;
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
				options->writeSiteSummary = true;
				break;
			case 'f':
			{
				options->filterAlignment = true;
				stringstream ss(optarg);

				while (!ss.eof())
				{
					char d = ss.get();
					switch(d)
					{
						case 'c':
							ss >> options->minCo;
							cout << "minCo=" << options->minCo << endl;
							break;
						case 'p':
							ss >> options->minPOC;
							cout << "minPOC=" << options->minPOC << endl;
							break;
						case 's':
							ss >> options->maxSmin;
							cout << "maxSmin=" << options->maxSmin << endl;
							break;
						case 'e':
							ss >> options->maxEntropy;
							cout << "maxEntropy=" << options->maxEntropy << endl;
							break;
						case ',':
							break;
						default:
							string param;
							while (ss.peek() != ',' && ss.peek() != EOF)
								param+= ss.get();
							cerr << "Unknown alignment filter " << d << param << "." << endl;
							break;
					}
				}
				break;
			}
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

	if (options->prefix.length() == 0)
	{
		int m = options->inputAlignment.find_last_of('/') + 1;
		int n = options->inputAlignment.find_last_of('.');
		if (n > -1)
			n = n-m;
		options->prefix = options->inputAlignment.substr(m, n);
	}

	return 0;
}

void printSyntax()
{
	cout << "Usage:" << endl;
	cout << "  shuffle [options] <alignment>" << endl;
	cout << "  shuffle -h" << endl;
	cout << endl;

	cout << "Options:" << endl;
	cout << "  -t<a|d|n>      Data type a=AA, d=DNA, n=Alphanumeric [default: auto-detect]" << endl;
	cout << "  -p<STRING>     Prefix for output files [default: name of alignment w/o .ext]" << endl;
	cout << "  -g<LIST>       Grouping of sites, e.g. 0,1,-2 for duplets, 0,1,2 for codons" << endl;
	cout << endl;
	cout << "  -d             Remove duplicates and write reduced alignment file" << endl;
	cout << "  -i             Remove informative sites duplicates, write reduced alignment" << endl;
	cout << endl;
	cout << "  -y             Perform tests of pairwise symmetry" << endl;
	cout << "  -x             Write extended output files with test results" << endl;
	cout << "  -w<NUM>,[NUM]  Window size and step width for symmetry test" << endl;
	cout << endl;
	cout << "  -r<NUM>        Number of randomizations for POC computations [default: 100]" << endl;
	cout << "  -s             Write a site summary" << endl;
	cout << "  -f<LIST>       Write a new alignment file, filtered by (comma-separated):" << endl;
	cout << "                   c<NUM>  Minimum Co score [default: 0]" << endl;
	cout << "                   p<NUM>  Minimum POC score [default: 0.0]" << endl;
	cout << "                   s<NUM>  Maximum Smin [default: " << INT_MAX << "]" << endl;
	cout << "                   e<NUM>  Maximum entropy [default: " << DBL_MAX << "]" << endl;
	cout << endl;
#ifdef _OPENMP
	cout << "  -n<NUM>        Number of threads [default: " << omp_get_max_threads() << "]" << endl;
#endif
	cout << "  -v[NUM]        Be increasingly verbose" << endl;
	cout << "  -h             This help page" << endl;
	cout << endl;

	cout << "Note:" << endl;
	cout << "  .phy or .phylip extensions indicate sequential Phylip file format" << endl;
	cout << "  .fsa or .fasta extensions indicate Fasta file format" << endl;


	exit(254);
}


int master(int argc, char** argv)
{
	Options options;

	cout << PROGNAME << " " << VERSION << "|";
#ifdef _OPENMP
	cout << "OpenMP|";
#elif _MPI
	cout << "MPI|";
#endif

	cout << PROGDATE << endl << endl;;

   	cout << "Parallel execution with " << numProcs << " processes." << endl << endl;

	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;

	if (!options.inputAlignment.length() || options.help)
		printSyntax();

#ifdef _MPI
	MPI_Bcast(&(options.randomizations), 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	Alignment alignment(&options);

	if (options.removeDuplicates)
	{
		alignment.removeDuplicates();
		alignment.write(options.prefix+".noDupes.fsa");
	}

	if (options.writeSiteSummary || options.filterAlignment || options.symmetryTest)
	{
		alignment.collectSites(&options);

		if (options.removeInformativeSitesDuplicates)
		{
			alignment.removeInformativeSitesDuplicates();
			alignment.write(options.prefix+".noDupes2.fsa");
		}

		if (options.symmetryTest)
			alignment.testSymmetry(options.prefix, options.writeExtendedTestResults, options.windowSize, options.windowStep);

		if (options.writeSiteSummary || options.filterAlignment)
		{
			alignment.computeBasicScores();
			alignment.send();
			alignment.computeCompatibilityScores(options.randomizations);
		}

		if (options.writeSiteSummary)
			alignment.writeSummary(options.prefix);

		if (options.filterAlignment)
		{
			cout << "Creating new alignment with minCo=" << options.minCo << " minPOC=" << options.minPOC << " maxSmin=" << options.maxSmin << " maxEntropy=" << options.maxEntropy << endl;
			Alignment modifiedAlignment = alignment.getModifiedAlignment(options.minCo, options.minPOC, options.maxSmin, options.maxEntropy);
			cout << "New alignment contains " << modifiedAlignment.getNumOfRows() << " sequences with " << modifiedAlignment.getNumOfCols() << " columns." << endl;
			modifiedAlignment.write(options.prefix+".filtered.fsa");
		}
	}

	return 0;
}


int worker()
{
	int randomizations;

	MPI_Bcast(&randomizations, 1, MPI_INT, 0, MPI_COMM_WORLD);

	Alignment alignment;
	alignment.recv();
	alignment.computeCompatibilityScores(randomizations);

	return 0;
}


int main(int argc, char** argv)
{
#ifdef _MPI
    if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
    {
    	cerr << "MPI_Init failed." << endl;
    	return -1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);
#elif _OPENMP
    numProcs = omp_get_max_threads();
#endif


    if (myId == 0)
    	master(argc, argv);
    else
    	worker();


#ifdef _MPI
	MPI_Finalize();
#endif

	return 0;
}
