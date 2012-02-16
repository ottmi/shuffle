#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <algorithm>
#include <climits>
#include <cfloat>
#include <string>
#include <sstream>
#include <cstdio>
#include "globals.h"
#include "Alignment.h"
#include "helper.h"
#ifdef _MPFR
#include "mpfr.h"
#endif

using namespace std;

int verbose = 0;
int numProcs = 1;

int parseArguments(int argc, char** argv, Options *options)
{
	char c;

	if (argc < 2)
		return 0;

	if (argc == 2 && argv[1][0] == '-')
	{
		options->help = true;
		return 0;
	}

	options->inputAlignment = string(argv[--argc]);

	options->alignmentFormat = -1;
	options->dataType = -1;
	options->columnFrom = -1;
	options->columnTo = -1;
	options->removeDuplicates = false;
	options->removeInformativeSitesDuplicates = false;
	options->writeInformativeSitesAlignment = false;
	options->removeIncompatibles = .0;
	options->convertAlignment = false;
	options->randomizations = 100;
	options->writeSiteSummary = false;
	options->writeRandomizedCo = false;
	options->filterAlignment = false;
	options->minCo = 0;
	options->minPOC = 0.0;
	options->maxSmin = INT_MAX;
	options->maxEntropy = DBL_MAX;
	options->help = 0;
	options->grouping.push_back(1);

	int minGroup = 0;
	int maxGroup = 0;

#ifdef _MPI
	string parameters = "t:p:c:g:deixr:sf:a:v::h";
#elif _OPENMP
	string parameters = "t:p:c:g:deij:xr:szf:a:n:v::h";
#else
	string parameters = "t:p:c:g:deij:xr:szf:a:v::h";
#endif

	while ((c = getopt(argc, argv, parameters.c_str())) != -1)
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
			case 'c':
			{
				stringstream ss(optarg);
				ss >> options->columnFrom;
				if (options->columnFrom == 0)
				{
					cerr << "Alignment column enumeration starts with 1" << endl;
					return 4;
				}
				if (ss.peek() == '-')
				{
					ss.ignore();
					ss >> options->columnTo;
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
					if (i == 0)
					{
						cerr << "The site numbers used in groupings start with 1." << endl;
						return 3;
					}
					if (i > 0)
						options->grouping.push_back(i);
					else
						i = -i;
					if (i > maxGroup) maxGroup = i;
					if (i < minGroup) minGroup = i;

					if (ss.peek() == ',') ss.ignore();
				}
				break;
			}
			case 'd':
				options->removeDuplicates = true;
				break;
			case 'e':
				options->removeInformativeSitesDuplicates = true;
				break;
			case 'i':
				options->writeInformativeSitesAlignment = true;
				break;
			case 'j':
			{
				stringstream ss(optarg);
				ss >> options->removeIncompatibles;
				break;
			}
			case 'x':
				options->convertAlignment = true;
				if (options->alignmentFormat == -1)
					options->alignmentFormat = -2;
				break;
			case 'r':
				options->randomizations = atoi(optarg);
				break;
			case 's':
				options->writeSiteSummary = true;
				break;
			case 'z':
				options->writeRandomizedCo = true;
				break;
			case 'f':
			{
				options->filterAlignment = true;
				stringstream ss(optarg);

				while (!ss.eof())
				{
					char d = ss.get();
					switch (d)
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
								param += ss.get();
							cerr << "Unknown alignment filter " << d << param << "." << endl;
							break;
					}
				}
				break;
			}
			case 'a':
			{
				char type = optarg[0];
				switch (type)
				{
					case 'f':
					case 'F':
						options->alignmentFormat = _FASTA_FORMAT;
						break;
					case 'p':
					case 'P':
						options->alignmentFormat = _PHYLIP_FORMAT;
						break;
					default:
						cerr << "Unknown alignment format: " << optarg << endl;
						return 3;
						break;
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

	if (options->columnFrom > 0) cout << "Columns: " << options->columnFrom << " - " << options->columnTo << endl;

	options->groupLength = maxGroup - minGroup + 1;

	if (options->prefix.length() == 0)
	{
		int m = options->inputAlignment.find_last_of('/') + 1;
		int n = options->inputAlignment.find_last_of('.');
		if (n > -1)
			n = n - m;
		options->prefix = options->inputAlignment.substr(m, n);
	}

	options->requireInformative = options->removeInformativeSitesDuplicates || options->writeInformativeSitesAlignment|| options->removeIncompatibles > 0 || options->writeSiteSummary || options->writeRandomizedCo || options->filterAlignment;

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
	cout << "  -c<from-to>    Only consider alignment columns from-to, enumeration starts with 1" << endl;
	cout << "  -g<LIST>       Grouping of sites, e.g. 0,1,-2 for duplets, 0,1,2 for codons" << endl;
	cout << endl;
	cout << "  -x             Convert alignment format" << endl;
	cout << "  -d             Remove duplicates and write reduced alignment file" << endl;
	cout << "  -e             Remove informative sites duplicates, write reduced alignment" << endl;
	cout << "  -i             Write reduced alignment only with parsimony informative sites" << endl;
#ifndef _MPI
	cout << "  -j<NUM>        Iteratively remove incompatible sites until avgCo>=NUM" << endl;
#endif
	cout << endl;
	cout << "  -r<NUM>        Number of randomizations for POC computations [default: 100]" << endl;
	cout << "  -s             Write a site summary" << endl;
#ifndef _MPI
	cout << "  -z             Write Co scores of randomized sites to file" << endl;
#endif
	cout << "  -f<LIST>       Write a new alignment file, filtered by (comma-separated):" << endl;
	cout << "                   c<NUM>  Minimum Co score [default: 0]" << endl;
	cout << "                   p<NUM>  Minimum POC score [default: 0.0]" << endl;
	cout << "                   s<NUM>  Maximum Smin [default: " << INT_MAX << "]" << endl;
	cout << "                   e<NUM>  Maximum entropy [default: " << DBL_MAX << "]" << endl;
	cout << endl;
	cout << "  -a<f|p>        Write [F]asta or [P]hylip alignments [default: same as input]" << endl;
#ifdef _OPENMP
	cout << "  -n<NUM>        Number of threads [default: " << omp_get_max_threads() << "]" << endl;
#endif
	cout << "  -v[NUM]        Be increasingly verbose" << endl;
	cout << "  -h             This help page" << endl;
	cout << endl;
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

#ifdef _MPFR
	cout << "MPFR|";
#endif
#ifdef _DEBUG
	cout << "DEBUG|";
#endif

	cout << PROGDATE << endl << endl;;

	int ret = parseArguments(argc, argv, &options);
	if (ret)
		return ret;

	if (!options.inputAlignment.length() || options.help)
	{
		printSyntax();
		return 254;
	}

#ifdef _MPFR
	cout << "Using MPFR " << MPFR_VERSION_STRING << " for high precision entropy computation." << endl;
#endif

#ifdef _MPI
   	cout << "Parallel execution with " << numProcs << " processes." << endl << endl;
	if (verbose)
	    cout << "This is master " << getMyId() << endl;
   	int buf[5];
   	buf[0] = verbose;
   	buf[1] = options.randomizations;
   	buf[2] = (int) options.writeRandomizedCo;
	MPI_Bcast(buf, 5, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	try
	{
		Alignment alignment = Alignment(&options);

		if (options.convertAlignment)
			alignment.write(options.prefix, options.alignmentFormat);

		if (options.removeDuplicates)
		{
			alignment.removeDuplicates();
			alignment.write(options.prefix + ".noDupes", options.alignmentFormat);
		}

		if (options.requireInformative)
		{
			alignment.collectSites(&options);

			if (verbose >= 4)
				alignment.printSites();

			if (options.writeInformativeSitesAlignment)
			{
			    Alignment informativeSitesAlignment = alignment.getInformativeSitesAlignment();
			    informativeSitesAlignment.write(options.prefix+".informative", options.alignmentFormat);
			}

			if (options.removeInformativeSitesDuplicates)
			{
				alignment.removeInformativeSitesDuplicates();
				alignment.write(options.prefix + ".noDupes2", options.alignmentFormat);
			}

			if (options.removeIncompatibles || options.writeSiteSummary || options.writeRandomizedCo || options.filterAlignment)
			{
				alignment.computeContextIndependentScores();
#ifdef _MPI
				alignment.send();
#endif
				alignment.computeContextDependentScores(options.randomizations, options.writeRandomizedCo);
			}

			if (options.writeSiteSummary)
				alignment.writeSummary(options.prefix + ".sites");

			if (options.filterAlignment)
			{
				Alignment filteredAlignment = alignment.getFilteredAlignment(options.minCo, options.minPOC, options.maxSmin, options.maxEntropy);
				cout << "New alignment contains " << filteredAlignment.getNumOfRows() << " sequences with " << filteredAlignment.getNumOfCols() << " columns." << endl;
				filteredAlignment.write(options.prefix + ".filtered", options.alignmentFormat);
			}
#ifndef _MPI
			if (options.writeRandomizedCo)
				alignment.writeRandomizedCo(options.prefix);

			if (options.removeIncompatibles > 0)
			    alignment.removeIncompatiblesIterative(&options);
#endif
		}

	} catch (string& s)
	{
		cerr << s << endl;
		return(255);
	}

	return 0;
}

#ifdef _MPI
int worker()
{
	int buf[5];

	MPI_Bcast(buf, 5, MPI_INT, 0, MPI_COMM_WORLD);
	verbose = buf[0];

	if (verbose)
	    cout << "This is worker " << getMyId() << endl;

	Alignment alignment;
	alignment.recv();
	alignment.computeContextDependentScores(buf[1], (bool) buf[2]);

	return 0;
}
#endif

int main(int argc, char** argv)
{

#ifdef _MPI
    if (MPI_Init(&argc,&argv)!=MPI_SUCCESS)
    {
    	cerr << "MPI_Init failed." << endl;
    	return -1;
    }

    numProcs = getNumOfCpus();
    if (getMyId() == 0)
    	master(argc, argv);
    else
    	worker();

    MPI_Finalize();
#else
    numProcs = getNumOfCpus();
    master(argc, argv);
#endif

    return 0;
}
