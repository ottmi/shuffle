#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include "FastaReader.h"
#include "PhylipReader.h"
#include "AASite.h"
#include "DNASite.h"
#include "AlphanumericSite.h"
#include "Alignment.h"

Alignment::Alignment(int dataType)
{
	_dataType = dataType;
}


Alignment::Alignment(Options *options)
{
	AlignmentReader *alignmentReader;

	string ext = options->inputAlignment.substr(options->inputAlignment.find_last_of('.') + 1);
	if (!ext.compare("phy") || !ext.compare("phylip"))
	{
		alignmentReader = new PhylipReader(options->inputAlignment);
	}
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
	{
		alignmentReader = new FastaReader(options->inputAlignment);
	}
	else
	{
		cerr << "Unknown input alignment format" << endl;
		exit(255);
	}
	_alignment = alignmentReader->getSequences();
	delete alignmentReader;

	string dataTypeDesc[] = {"DNA", "AA", "alphanumeric"};
	if (options->dataType < 0)
	{
		map<char,unsigned long> baseOccurences;
		for (unsigned int i=0; i<_alignment.size(); i++)
		{
			string s = _alignment[i].getSequence();
			for (unsigned int j=0; j<s.length(); j++)
				baseOccurences[s[j]]++;
		}

		string maps[] = {_DNA_MAP, _AA_MAP, _ALPHANUM_MAP};
		unsigned long counts[3];
		for (unsigned int i=0; i<3; i++)
		{
			counts[i] = 0;
			string map = maps[i];
			for (unsigned j=0; j<map.length(); j++)
				counts[i]+= baseOccurences[map[j]];
		}

		if (verbose)
			cout << counts[0] << " DNA characters, " << counts[1] << " AA characters, " << counts[2] << " Alphanum characters." << endl;
		int dataTypeGuess = _ALPHANUM_DATA;
		if (counts[2] == counts[0])
			dataTypeGuess = _DNA_DATA;
		else if (counts[2] == counts[1])
			dataTypeGuess = _AA_DATA;
		_dataType = dataTypeGuess;

		cout << "The alignment contains " << getNumOfRows() << " sequences " << "which appear to be " << dataTypeDesc[_dataType] << "." << endl;
	} else
	{
		_dataType = options->dataType;
		cout << "The alignment contains " << getNumOfRows() << " sequences which have been defined to be " << dataTypeDesc[_dataType] << "." << endl;
	}
}


Alignment::~Alignment()
{
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
		delete _informativeSites[i];
}


void Alignment::addSequence(Sequence s)
{
	_alignment.push_back(s);
}


void Alignment::removeDuplicates()
{
	cout << endl;
	cout << "Removing duplicates...";
	if (verbose)
		cout << endl;
	vector<Sequence>::iterator it1, it2;
	int count = 0;
	for (it1=_alignment.begin(); it1!=_alignment.end(); it1++)
	{
		it2 = it1+1;
		while (it2 != _alignment.end())
		{
			if (it1->getSequence() == it2->getSequence())
			{
				if (verbose)
					cout << "  " << it2->getName() << " is a duplicate of " << it1->getName() << endl;
				_alignment.erase(it2);
				count++;
			}	else
			{
				it2++;
			}
		}
	}
	if (!verbose)
		cout << "\b\b\b, done." << endl;

	cout << "Removed " << count << " duplicates, " << getNumOfRows() << " sequences remain in the alignment." << endl;
}


void Alignment::collectInformativeSites(Options *options)
{
	cout << endl;
	cout << "Collecting informative sites..." << endl;
	long t1 = time(NULL);
	unsigned int numOfSites = (getNumOfCols()-options->groupOffset) / options->groupLength;

#ifdef _OPENMP
	_informativeSites.resize(numOfSites, NULL);
	#pragma omp parallel for
#endif
	for (unsigned int i = 0; i < numOfSites; i++)
	{
		Site *s = NULL;
		switch (_dataType)
		{
			case _DNA_DATA:
				s = new DNASite(&_alignment, options->grouping, options->groupLength*i);
				break;
			case _AA_DATA:
				s = new AASite(&_alignment, options->grouping, options->groupLength*i);
				break;
			case _ALPHANUM_DATA:
				s = new AlphanumericSite(&_alignment, options->grouping, options->groupLength*i);
				break;
			default:
				cerr << "Unknown data type " << _dataType << " at Alignment::Alignment()" << endl;
				break;
		}

		if (s)
		{
			if (s->isInformative())
			{
				#ifdef _OPENMP
					_informativeSites[i] = s;
				#else
					_informativeSites.push_back(s);
				#endif
			}	else
			{
					delete s;
			}
		}
	}

#ifdef _OPENMP
	vector<Site*>::iterator it=_informativeSites.begin();
	while (it != _informativeSites.end())
	{
		if (*it == NULL)
			_informativeSites.erase(it);
		else
			it++;
	}
#endif
	long t2 = time(NULL);

	cout << "Found " << numOfSites << " sites, " << _informativeSites.size() << " of which are informative, taking " << t2-t1 << "s" << endl;

	if (options->groupLength > 1)
		cout << "Each site consists of " << options->groupLength << " out of " << getNumOfCols() << " columns, starting with an offset of " << options->groupOffset << "." << endl;
}

void Alignment::computeCompatibilityScores(int randomizations)
{
	cout << endl;
	cout << "Computing compatibility scores, doing " << randomizations << " randomizations..." << endl;
	cout << "0%" << flush;

	long t1 = time(NULL);
	unsigned long n = _informativeSites.size();
	unsigned long total = n*(n-1)/2 + n*n*randomizations;
	unsigned long count = 0;
	int myTid = 0;

#ifdef _OPENMP
#pragma omp parallel shared(count) private(myTid)
	{
		myTid = omp_get_thread_num();
#pragma omp for
#endif
	for (unsigned int i = 0; i < n; i++)
	{
		for (unsigned int j = i + 1; j < n; j++)
		{
			if (_informativeSites[i]->checkCompatibility(_informativeSites[j]))
			{
				_informativeSites[i]->incComp();
				_informativeSites[j]->incComp();
			}
		}
		count+= n-i-1;
		if (myTid == 0)
		{
			long elapsed = time(NULL) - t1;
			cout << "\r" << count * 100 / total << "%\tTime elapsed: " << elapsed/60 << ":" << setfill('0') << setw(2) << elapsed%60 << flush;
		}
	}

#ifdef _OPENMP
#pragma omp for
#endif
	for (unsigned int i = 0; i < n; i++)
	{
		_informativeSites[i]->computeScores(n);

		if (randomizations)
		{
			int poc = 0;
			for (int r = 0; r < randomizations; r++)
			{
				int comp = 0;
				Site* randomSite = _informativeSites[i]->randomize();
				for (unsigned int j = 0; j < n; j++)
				{
					if (i != j && randomSite->checkCompatibility(_informativeSites[j]))
						comp++;
				}
				delete randomSite;
				if (_informativeSites[i]->getComp() <= comp)
					poc++;
			}

			count += randomizations*n;
			if (myTid == 0)
			{
				long elapsed = time(NULL) - t1;
				long eta = (elapsed * total) / count - elapsed;
				cout << "\r" << count * 100 / total << "%\tTime elapsed: " << elapsed/60 << ":" << setfill('0') << setw(2) << elapsed%60 << "\tETA: " << eta/60 << ":" << setw(2) << eta%60 << "  " << flush;
			}
			_informativeSites[i]->setPOC((double) poc / randomizations);
		}
		else
			_informativeSites[i]->setPOC(.0);
	}
#ifdef _OPENMP
}
#endif
	long t2 = time(NULL);
	cout << "\rFinished computing scores, taking " << t2 - t1 << "s.               " << endl;
}


Alignment Alignment::getModifiedAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy)
{
	Alignment a(_dataType);
	cout << "Creating new alignment with minCo=" << minCo << " minPOC=" << minPOC << " maxSmin=" << maxSmin << " maxEntropy=" << maxEntropy << endl;
	vector<int> sites;

	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site *site = _informativeSites[i];
		if (site->getCo() >= minCo && site->getPOC() >= minPOC && site->getSmin() <= maxSmin && site->getEntropy() <= maxEntropy)
			sites.push_back(i);
	}

	for (unsigned int i = 0; i < _alignment.size(); i++)
	{
		string newSeq;
		Sequence seq = _alignment[i];
		for (unsigned int j = 0; j < sites.size(); j++)
		{
			Site *site = _informativeSites[sites[j]];
			newSeq += seq.getColumns(site->getCols());
		}
		if (newSeq.length())
		{
			Sequence s(_alignment[i].getName(), newSeq);
			a.addSequence(s);
		}
	}
	cout << "New alignment contains " << a.getNumOfRows() << " sequences with " << a.getNumOfCols() << " columns." << endl;

	return a;
}


void Alignment::writeSummary(string fileName)
{
	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("\n\nError, cannot open file " + fileName);
	cout << "Writing site summary to " << fileName << endl;

	file << "Site No.,Smin,Entropy,OV,Co,Poc" << endl;
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		file << s->getCols()[0] + 1 << "," << s->getSmin() << "," << s->getEntropy() << "," << s->getOV() << "," << s->getCo() << "," << s->getPOC() << endl;
	}
}


void Alignment::write(string fileName)
{
	int type;
	string ext = fileName.substr(fileName.find_last_of('.') + 1);
	if (!ext.compare("phy") || !ext.compare("phylip"))
		type = 0;
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
		type = 1;
	else
	{
		cerr << "Couldn't recognize output format." << endl;
		cerr << "Please use extension .phy or .phylip for Phylip format," << endl;
		cerr << "and .fsa or .fasta for Fasta format." << endl;
		return;
	}

	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("\n\nError, cannot open file " + fileName);

	switch (type)
	{
		case 0:
			cout << "Writing Phylip alignment to " << fileName << endl;
			file << getNumOfRows() << " " << getNumOfCols() << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
				file << _alignment[i].getName() << " " << _alignment[i].getSequence() << endl;
			break;

		case 1:
			cout << "Writing Fasta alignment to " << fileName << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
			{
				file << ">" << _alignment[i].getName() << endl;
				file << _alignment[i].getSequence() << endl;
			}
			break;
	}
	file.close();
}
