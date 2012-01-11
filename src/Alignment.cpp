#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <set>
#include <cmath>
#include <climits>
#include <cfloat>
#include <sstream>
#include "AlignmentReader.h"
#include "AASite.h"
#include "DNASite.h"
#include "AlphanumericSite.h"
#include "Alignment.h"
#include "helper.h"


Alignment::Alignment()
{
	_cols = 0;
}


Alignment::Alignment(Options *options)
{
	AlignmentReader alignmentReader(options->inputAlignment);
 	_alignment = alignmentReader.getSequences();
 	_cols = alignmentReader.getCols();

 	if (options->alignmentFormat == -1)
 		options->alignmentFormat = alignmentReader.getFormat();
 	else if (options->alignmentFormat == -2)
 		options->alignmentFormat = (alignmentReader.getFormat() + 1)%2;

	string dataTypeDesc[] = { "DNA", "AA", "alphanumeric" };
	if (options->dataType < 0)
	{
		map<char, unsigned long> baseOccurences;
		for (unsigned int i = 0; i < _alignment.size(); i++)
		{
			string s = _alignment[i].getSequence();
			for (unsigned int j = 0; j < s.length(); j++)
				baseOccurences[s[j]]++;
		}

		string maps[] = { _DNA_MAP, _AA_MAP, _ALPHANUM_MAP };
		unsigned long counts[3];
		for (unsigned int i = 0; i < 3; i++)
		{
			counts[i] = 0;
			string map = maps[i];
			for (unsigned j = 0; j < map.length(); j++)
				counts[i] += baseOccurences[map[j]];
		}

		if (verbose)
			cout << counts[0] << " DNA characters, " << counts[1] << " AA characters, " << counts[2] << " Alphanum characters." << endl;
		int dataTypeGuess = _ALPHANUM_DATA;
		if (counts[2] == counts[0])
			dataTypeGuess = _DNA_DATA;
		else if (counts[2] == counts[1])
			dataTypeGuess = _AA_DATA;
		_dataType = dataTypeGuess;

		cout << "It contains " << getNumOfRows() << " sequences " << "which appear to be " << dataTypeDesc[_dataType] << "." << endl;
	} else
	{
		_dataType = options->dataType;
		cout << "It contains " << getNumOfRows() << " sequences which have been defined to be " << dataTypeDesc[_dataType] << "." << endl;
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
	if (s.getLength() > _cols)
		_cols = s.getLength();
}


void Alignment::removeDuplicates()
{
	cout << endl;
	cout << "Removing duplicates...";
	if (verbose)
		cout << endl;
	vector<Sequence>::iterator it1, it2;
	int count = 0;
	for (it1 = _alignment.begin(); it1 != _alignment.end(); it1++)
	{
		it2 = it1 + 1;
		while (it2 != _alignment.end())
		{
			if (it1->getSequence() == it2->getSequence())
			{
				if (verbose)
					cout << "  " << it2->getName() << " is a duplicate of " << it1->getName() << endl;
				it2 = _alignment.erase(it2);
				count++;
			} else
			{
				it2++;
			}
		}
	}
	if (!verbose)
		cout << "\b\b\b, done." << endl;

	cout << "Removed " << count << " duplicates, " << getNumOfRows() << " sequences remain in the alignment." << endl;
}


void Alignment::removeInformativeSitesDuplicates()
{
	cout << endl;
	cout << "Removing duplicates based only on informative sites...";
	if (verbose)
		cout << endl;

	vector<unsigned int> eraseList;
	unsigned int total = 0;
	unsigned int count = 0;
	do
	{
		vector<Site*> sites = _informativeSites;
		vector<Sequence> a = getSubAlignment(sites).getAlignment();

		vector<Sequence>::iterator it1, it2;
		unsigned int i = 0;
		count = 0;

		for (it1 = a.begin(); it1 != a.end(); it1++)
		{
			it2 = it1 + 1;
			unsigned j = i + 1;
			while (it2 != a.end())
			{
				if (it1->getSequence() == it2->getSequence())
				{
					if (verbose)
						cout << "  " << it2->getName() << " is a duplicate of " << it1->getName() << endl;
					count++;
					eraseList.push_back(j);
					it2 = a.erase(it2);
				} else
				{
					it2++;
				}
				j++;
			}

			if (eraseList.size() > 0)
			{
				for (unsigned int k = eraseList.size(); k--; )
				{
					unsigned int l = eraseList[k];
					for (unsigned int m = 0; m < _sites.size(); m++)
						_sites[m]->remove(l);
					_alignment.erase(_alignment.begin()+l);
				}
				eraseList.clear();
			}
			i++;
		}

		_informativeSites.clear();
		for (unsigned int m = 0; m < _sites.size(); m++)
		{
			if (_sites[m]->checkInformative())
				_informativeSites.push_back(_sites[m]);
		}
		total+= count;
	} while (count);

	if (!verbose)
		cout << "\b\b\b, done." << endl;

	cout << "Removed " << total << " duplicates, " << getNumOfRows() << " sequences with " << _informativeSites.size() << " informative sites remain in the alignment." << endl;
}

#ifndef _MPI
bool compPos(Site *a, Site *b)
{
	return a->getCols()[0] < b->getCols()[0];
}

bool compCo(Site *a, Site *b)
{
	unsigned int aVal = a->getComp();
	unsigned int bVal = b->getComp();

	if (aVal == bVal)
		return compPos(a, b);
	else
		return aVal < bVal;
}

void Alignment::removeIncompatiblesIterative(Options *options)
{
	double threshold = options->removeIncompatibles;
	cout << endl << "Removing incompatbile sites iteratively, target avgCo=" << threshold << endl;

	double avgCo = .0;
	double siteCo;
	unsigned int siteComp;
	int siteCol;
	while (avgCo < threshold)
	{
		sort(_informativeSites.begin(), _informativeSites.end(), compCo);
		Site *s = _informativeSites.front();
		siteCol = s->getCols()[0];
		siteCo = s->getCo();
		siteComp = s->getComp();

		_informativeSites.erase(_informativeSites.begin());

		avgCo = .0;
		for (unsigned int i = 0; i < _informativeSites.size(); i++)
		{
			_informativeSites[i]->removeCompatibleSite(siteCol);
			_informativeSites[i]->computeCo(_informativeSites.size());
			avgCo+= _informativeSites[i]->getCo();
		}
		avgCo /= _informativeSites.size();

		if (verbose)
			cout << "  Site " << siteCol+1 << ": Comp=" << siteComp << " Co=" << siteCo << " avgCo=" << avgCo << endl;
	}

	sort(_informativeSites.begin(), _informativeSites.end(), compPos);

	stringstream ss;
	ss << options->prefix << ".iterative." << threshold;
	writeSummary(ss.str());
	Alignment a = getSubAlignment(_informativeSites);
	a.write(ss.str(), options->alignmentFormat);
}
#endif

void Alignment::collectSites(Options *options)
{
	cout << endl;
	cout << "Collecting sites..." << endl;
	long t1 = time(NULL);
	long lastTime = t1;
	unsigned int numOfSites = (getNumOfCols() - options->groupOffset) / options->groupLength;
	unsigned int count = 0;
	_sites.resize(numOfSites, NULL);

#ifdef _OPENMP
#pragma omp parallel for shared (count) schedule(guided)
#endif
	for (unsigned int i = 0; i < numOfSites; i++)
	{
		Site *s = NULL;
		switch (_dataType) {
			case _DNA_DATA:
				s = new DNASite(i, options);
				break;
			case _AA_DATA:
				s = new AASite(i, options);
				break;
			case _ALPHANUM_DATA:
				s = new AlphanumericSite(i, options);
				break;
			default:
				cerr << "Unknown data type " << _dataType << " at Alignment::Alignment()" << endl;
				break;
		}
		if (s)
		{
			_sites[i] = s;
			s->initialize(&_alignment);
			if (options->requireInformative)
				s->checkInformative();
		}
		count++;
		if (myId == 0)
		{
			long t2 = time(NULL);
			if (t2 > lastTime)
			{
				long elapsed = t2 - t1;
				long eta = (elapsed * numOfSites) / count - elapsed;
				lastTime = t2;
				cout << "\r" << count * 100 / numOfSites << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;
			}
		}
	}

	long t2 = time(NULL);
	cout << "\rDone, taking " << printTime(t2-t1) << "                         " << endl;

	if (options->requireInformative)
	{
		for (unsigned int i = 0; i < _sites.size(); i++)
		{
			Site *s = _sites[i];
			if (s->isInformative())
				_informativeSites.push_back(s);
		}
	}

	if (options->requireInformative)
		cout << "Found " << numOfSites << " sites, " << _informativeSites.size() << " of which are informative." << endl;
	else
		cout << "Found " << numOfSites << " sites." << endl;

	if (options->groupLength > 1)
		cout << "Each site consists of " << options->groupLength << " out of " << getNumOfCols() << " columns, starting with an offset of " << options->groupOffset << "." << endl;
}


void Alignment::checkIdenticalSites()
{
	cout << endl;
	cout << "Checking for identical sites..." << endl;

	unsigned int n = _informativeSites.size();
	unsigned int count = 0;
	vector<bool> duplicates (n, false);
#ifdef _OPENMP
#pragma omp parallel for shared(count) schedule(guided)
#endif
	for (unsigned int i = 0; i < n; i++)
		for (unsigned int j = i + 1; j < n; j++)
			if (!duplicates[j] && _informativeSites[i]->compare(_informativeSites[j]))
			{
				cout << "  " << _informativeSites[i]->colsToString() << " and " << _informativeSites[j]->colsToString() << " are identical" << endl;
				duplicates[j] = true;
				count++;
			}

	cout << "Done, found " << count << " identical sites." << endl;
}


void Alignment::computeContextIndependentScores()
{
	cout << endl;
	cout << "Computing context-independent scores..." << endl;

	unsigned long n = _informativeSites.size();
	long t1 = time(NULL);


#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (unsigned int i = 0; i < n; i++)
		_informativeSites[i]->computeScores(n);

	long t2 = time(NULL);
	cout << "\rDone, taking " << printTime(t2-t1) << "                         " << endl;
}


void Alignment::computeCo(unsigned int start, unsigned int stop, unsigned int n)
{
	if (myId == 0)
	    cout << "  Computing Co:  0%" << flush;

	long count, total, t1, t2, lastTime;
	count = 0;
	total = stop - start;
	t1 = time(NULL);
	lastTime = t1;

#ifdef _OPENMP
	long chunk = n / (omp_get_num_threads() * 8);
#pragma omp parallel for shared(count) schedule(dynamic, chunk)
#endif
	for (unsigned int i = start; i <= stop; i++)
	{
		for (unsigned int j = i + 1; j < n; j++)
		{
			if (_informativeSites[i]->checkCompatibility(_informativeSites[j]))
			{
				_informativeSites[i]->addCompatibleSite(_informativeSites[j]->getCols()[0]);
				_informativeSites[j]->addCompatibleSite(_informativeSites[i]->getCols()[0]);
			}
		}

		count++;
		if (myId == 0)
		{
			t2 = time(NULL);
			if (t2 > lastTime)
			{
				long elapsed = t2 - t1;
				long eta = (elapsed * total) / count - elapsed;
				cout << "\r  Computing Co:  " << count * 100 / total << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;

			}
		}
	}

#ifndef _MPI
#ifdef _OPENMP
#pragma omp parallel for shared(count) schedule(guided)
#endif
	for (unsigned int i = start; i <= stop; i++)
		_informativeSites[i]->computeCo(n);
#endif
	if (myId == 0)
	{
	    t2 = time(NULL);
	    cout << "\r  Computing Co:  Done, taking " << printTime(t2-t1) << "                         " << endl;
	}
}


void Alignment::computePOC(unsigned int start, unsigned int stop, unsigned int n, unsigned int randomizations)
{
	if (myId == 0)
	    cout << "  Computing POC: 0%" << flush;

	long count, total, t1, t2, lastTime;
	count = 0;
	total = stop - start;
	t1 = time(NULL);
	lastTime = t1;

#ifdef _OPENMP
#pragma omp parallel for shared(count) schedule(guided)
#endif
	for (unsigned int i = start; i <= stop; i++)
	{
#ifdef _DEBUG
		srand( i+42 );
#endif
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
			_informativeSites[i]->addRandomizedCo(((double) comp) / n);
			if (_informativeSites[i]->getComp() <= comp)
				poc++;
		}

		count++;
		if (myId == 0)
		{
			t2 = time(NULL);
			if (t2 > lastTime)
			{
				long elapsed = t2 - t1;
				long eta = (elapsed * total) / count - elapsed;
				cout << "\r  Computing POC: " << count * 100 / total << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;
			}
		}
		_informativeSites[i]->computePOC(poc, randomizations);
	}

	if (myId == 0)
	{
	    t2 = time(NULL);
	    cout << "\r  Computing POC: Done, taking " << printTime(t2-t1) << "                         " << endl;
	}
}


void Alignment::computeR(unsigned int start, unsigned int stop, unsigned int n)
{
	if (myId == 0)
	    cout << "  Computing r: 0%" << flush;

	long count, total, t1, t2, lastTime;
	count = 0;
	total = stop - start;
	t1 = time(NULL);
	lastTime = t1;

#ifdef _OPENMP
#pragma omp parallel for shared(count) schedule(guided)
#endif
	for (unsigned int i = start; i <= stop; i++)
	{
	    double sum = 0;
	    for (unsigned int j = 0; j < n; j++)
	    {
		if (i != j)
		    sum+= _informativeSites[i]->checkPattern(_informativeSites[j]);
	    }
	    _informativeSites[i]->setR(sum/(n-1));

	    count++;
	    if (myId == 0)
	    {
		t2 = time(NULL);
		if (t2 > lastTime)
		{
		    long elapsed = t2 - t1;
		    long eta = (elapsed * total) / count - elapsed;
		    cout << "\r  Computing r: " << count * 100 / total << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;
		}
	    }
	}

	if (myId == 0)
	{
	    t2 = time(NULL);
	    cout << "\r  Computing r: Done, taking " << printTime(t2-t1) << "                         " << endl << endl;
	}
}

void Alignment::computeContextDependentScores(int randomizations)
{
	if (myId == 0)
	{
	    cout << endl;
	    cout << "Computing context-dependent scores, doing " << randomizations << " randomizations for POC:" << endl;
	}

	unsigned int start, end;
	unsigned int n = _informativeSites.size();
	srand( time(NULL) );
#ifdef _MPI
	unsigned int share = (n + numProcs - 1)/ numProcs;
	start = myId * share;
	end = (myId + 1) * share - 1;
	if (end >= _informativeSites.size())
		end = _informativeSites.size() - 1;

	computeCo(start, end, n);
	int *sendBuf1 = (int *) malloc(sizeof(int) * n);
	int *recvBuf1 = (int *) malloc(sizeof(int) * n);
	for (unsigned int i = 0; i < n; i++)
		sendBuf1[i] = _informativeSites[i]->getComp();
	MPI_Allreduce(sendBuf1, recvBuf1, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	for (unsigned int i = 0; i < n; i++)
	{
		_informativeSites[i]->setComp(recvBuf1[i]);
		_informativeSites[i]->computeCo(n);
	}
	free(recvBuf1);
	free(sendBuf1);

	if (randomizations)
	    computePOC(start, end, n, randomizations);

	for (unsigned int i = 0; i < n; i++)
	    _informativeSites[i]->checkInformative();
	computeR(start, end, n);

	unsigned int k=0;
	double *sendBuf2 = (double *) malloc(sizeof(double) * share * 2);
	for (unsigned int i = start; i <= end; i++)
	{
	    sendBuf2[k*2] = _informativeSites[i]->getPOC();
	    sendBuf2[k*2 + 1] = _informativeSites[i]->getR();
	    k++;
	}

	double *recvBuf2 = NULL;
	if (myId == 0)
		recvBuf2 = (double *) malloc(sizeof(double) * share * 2 * numProcs);
	MPI_Gather(sendBuf2, share*2, MPI_DOUBLE, recvBuf2, share*2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myId == 0)
	{
		for (unsigned int i = 0; i < n; i++)
		{
			_informativeSites[i]->setPOC(recvBuf2[i*2]);
			_informativeSites[i]->setR(recvBuf2[i*2 + 1]);
		}
		free(recvBuf2);
	}

	free(sendBuf2);
#else
	start = 0;
	end = n-1;

	computeCo(start, end, n);

	if (randomizations)
	    computePOC(start, end, n, randomizations);

	computeR(start, end, n);
#endif
}


Alignment Alignment::getFilteredAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy)
{
	cout << "Creating new alignment with minCo=" << minCo << " minPOC=" << minPOC << " maxSmin=" << maxSmin << " maxEntropy=" << maxEntropy << endl;

	vector<Site*> sites;
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site *site = _informativeSites[i];
		if (site->getCo() >= minCo && site->getPOC() >= minPOC && site->getSmin() <= maxSmin && (isnan(site->getEntropy()) || site->getEntropy() <= maxEntropy))
			sites.push_back(site);
	}

	return getSubAlignment(sites);
}


Alignment Alignment::getInformativeSitesAlignment()
{
    return getSubAlignment(_informativeSites);
}


Alignment Alignment::getSubAlignment(vector<Site*> sites)
{
	Alignment a;
	for (unsigned int i = 0; i < _alignment.size(); i++)
	{
		string newSeq;
		Sequence seq = _alignment[i];
		vector<Site*>::iterator it;
		for (it = sites.begin(); it != sites.end(); it++)
			newSeq += seq.getColumns((*it)->getCols());

		if (newSeq.length())
		{
			Sequence s(_alignment[i].getName(), newSeq);
			a.addSequence(s);
		}
	}

	return a;
}


void Alignment::writeRandomizedCo(string prefix)
{
	string fileName = prefix + ".poc.csv";
	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("Error, cannot open file " + fileName);
	cout << "Writing Co scores of randomized sites to " << fileName << endl;

	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		vector<double> co = s->getRandomizedCo();
		file << s->getCols()[0] + 1;
		for (unsigned int j = 0; j < co.size(); j++)
			file << "," << scientific << co[j];
		file << endl;
	}
	file.close();
}


void Alignment::writeSummary(string prefix)
{
	string fileName = prefix + ".csv";
	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("Error, cannot open file " + fileName);
	cout << "Writing site summary to " << fileName << endl;

	set<int> bases;
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		BaseOccurenceMap f = s->getFrequencies();
		for (BaseOccurenceMapIterator it = f.begin(); it != f.end(); it++)
			bases.insert(it->first);
	}

	file << "Site No.,Smin,Entropy,OV,Co,Poc,r";
	Site* s = _informativeSites[0];
	for (set<int>::iterator it = bases.begin(); it != bases.end(); it++)
		file << ",f(" << s->mapNumToChar(*it) << ")";
	file << ",f(?)" << endl;

	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		BaseOccurenceMap f = s->getFrequencies();
		file << s->getCols()[0] + 1 << "," << s->getSmin() << "," << fixed << s->getEntropy() << "," << s->getOV() << "," << s->getCo() << "," << s->getPOC()<< "," << s->getR();
		for (set<int>::iterator it = bases.begin(); it != bases.end(); it++)
			file << "," << (((double) f[*it]) / getNumOfRows());
		file << "," << ((double) s->getAmbiguousCount()) / getNumOfRows() << endl;
	}
}


void Alignment::write(string baseName, int format)
{
	string fileName;
	if (format == 0)
		fileName = baseName + ".fsa";
	else
		fileName = baseName + ".phy";

	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("Error, cannot open file " + fileName);

	switch (format) {
		case _FASTA_FORMAT:
			cout << "Writing Fasta alignment to " << fileName << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
			{
				file << ">" << _alignment[i].getName() << endl;
				file << _alignment[i].getSequence() << endl;
			}
			break;

		case _PHYLIP_FORMAT:
			cout << "Writing Phylip alignment to " << fileName << endl;
			file << getNumOfRows() << " " << getNumOfCols() << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
			{
				string name = _alignment[i].getName();
				if (name.length() >= 10)
					file << name << " ";
				else
					file << setw(10) << left << name;
				file << _alignment[i].getSequence() << endl;
			}
			break;
	}
	file.close();
}

#ifdef _MPI
void Alignment::send()
{
	unsigned int m = _informativeSites.size();
	unsigned int n = getNumOfRows();
	int buf[3];
	buf[0] = _dataType;
	buf[1] = (int) m;
	buf[2] = (int) n;

	unsigned int *buf2 = (unsigned int *) malloc(sizeof(unsigned int) * m * n);
	unsigned int *buf3 = (unsigned int *) malloc(sizeof(unsigned int) * m);
	for (unsigned int i = 0; i < m; i++)
	{
		Site *site = _informativeSites[i];
		vector<unsigned int> v = site->getSite();
		copy(v.begin(), v.end(), buf2+i*n);
		buf3[i] = site->getCols()[0];
	}

	MPI_Bcast(buf, 3, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(buf2, m*n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(buf3, m, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	free(buf2);
	free(buf3);
}

void Alignment::recv()
{
	int buf[3];
	int m, n;

	MPI_Bcast(buf, 3, MPI_INT, 0, MPI_COMM_WORLD);
	_dataType = buf[0];
	m = buf[1];
	n = buf[2];

	unsigned int *buf2 = (unsigned int *) malloc(sizeof(unsigned int) * m * n);
	MPI_Bcast(buf2, m*n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	unsigned int *buf3 = (unsigned int *) malloc(sizeof(unsigned int) * m);
	MPI_Bcast(buf3, m, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	vector<unsigned int> v;
	v.resize(n);
	for (unsigned int i = 0; i < m; i++)
	{
		copy(buf2+i*n, buf2+(i+1)*n, v.begin());
		Site *site;
		switch (_dataType)
		{
			case _DNA_DATA:
				site = new DNASite(buf3[i], v);
				break;
			case _AA_DATA:
				site = new AASite(buf3[i], v);
				break;
			case _ALPHANUM_DATA:
				site = new AlphanumericSite(buf3[i], v);
				break;
		}
		_informativeSites.push_back(site);
	}

	free(buf2);
	free(buf3);
}
#endif
