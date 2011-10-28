#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>
#include <algorithm>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "AlignmentReader.h"
#include "AASite.h"
#include "DNASite.h"
#include "AlphanumericSite.h"
#include "Alignment.h"
#include "helper.h"

Alignment::Alignment(int dataType)
{
	_dataType = dataType;
}

Alignment::Alignment(Options *options)
{
	AlignmentReader alignmentReader(options->inputAlignment);
 	_alignment = alignmentReader.getSequences();

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
		vector<Sequence> a = getModifiedAlignment(0, 0, INT_MAX, DBL_MAX).getAlignment();
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


void Alignment::collectSites(Options *options)
{
	cout << endl;
	cout << "Collecting sites..." << endl;
	long t1 = time(NULL);
	long lastTime = t1;
	unsigned int numOfSites = (getNumOfCols() - options->groupOffset) / options->groupLength;
	unsigned int count = 0;
	bool requireInformative = options->writeSiteSummary || options->filterAlignment;
	_sites.resize(numOfSites, NULL);

#ifdef _OPENMP
#pragma omp parallel for shared (count)
#endif
	for (unsigned int i = 0; i < numOfSites; i++)
	{
		Site *s = NULL;
		switch (_dataType) {
			case _DNA_DATA:
				s = new DNASite(&_alignment, i, options);
				break;
			case _AA_DATA:
				s = new AASite(&_alignment, i, options);
				break;
			case _ALPHANUM_DATA:
				s = new AlphanumericSite(&_alignment, i, options);
				break;
			default:
				cerr << "Unknown data type " << _dataType << " at Alignment::Alignment()" << endl;
				break;
		}
		if (s)
		{
			_sites[i] = s;
			if (requireInformative)
				s->checkInformative();
		}
		count++;
		if (omp_get_thread_num() == 0)
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

	if (requireInformative)
	{
		for (unsigned int i = 0; i < _sites.size(); i++)
		{
			Site *s = _sites[i];
			if (s->isInformative())
				_informativeSites.push_back(s);
		}
	}

	if (requireInformative)
		cout << "Found " << numOfSites << " sites, " << _informativeSites.size() << " of which are informative." << endl;
	else
		cout << "Found " << numOfSites << " sites." << endl;

	if (options->groupLength > 1)
		cout << "Each site consists of " << options->groupLength << " out of " << getNumOfCols() << " columns, starting with an offset of " << options->groupOffset << "." << endl;
}

void Alignment::testSymmetry(string prefix, bool extended, int windowSize, int windowStep)
{
	string resultsFileName = prefix + ".symmetry.csv";
	string bowkerFileName = prefix + ".bowker.csv";
	string delta_sFileName = prefix + ".delta_s.csv";
	string delta_msFileName = prefix + ".delta_ms.csv";
	ofstream resultsFile, bowkerFile, delta_sFile, delta_msFile;
	cout << endl << "Performing tests of pairwise symmetry, writing results to: ";
	if (extended)
	{
		cout << endl;
		cout << "  comprehensive spreadsheet: " << resultsFileName << endl;
		cout << "  Bowker matrix:             " << bowkerFileName << endl;
		cout << "  delta_s distance matrix:   " << delta_sFileName << endl;
		cout << "  delta_ms distance matrix:  " << delta_msFileName << endl;

		bowkerFile.open(bowkerFileName.c_str(), ifstream::trunc);
		if (!bowkerFile.is_open())
			throw("\n\nError, cannot open file " + bowkerFileName);

		delta_sFile.open(delta_sFileName.c_str(), ifstream::trunc);
		if (!delta_sFile.is_open())
			throw("\n\nError, cannot open file " + delta_sFileName);

		delta_msFile.open(delta_msFileName.c_str(), ifstream::trunc);
		if (!delta_msFile.is_open())
			throw("\n\nError, cannot open file " + delta_msFileName);
	} else
	{
		cout << resultsFileName << endl;
	}

	resultsFile.open(resultsFileName.c_str(), ifstream::trunc);
	if (!resultsFile.is_open())
		throw("\n\nError, cannot open file " + resultsFileName);

	unsigned int n = _alignment.size();
	unsigned int cols = _sites.size();
	if (windowSize <= 0)
		windowSize = cols;
	if (windowStep <= 0)
		windowStep = windowSize;

	if (windowSize < (int) cols)
		cout << "  WindowSize=" << windowSize << " StepWidth=" << windowStep << endl;

	int dim;
	if (_dataType == _DNA_DATA)
		dim = 5;
	else if (_dataType == _AA_DATA)
		dim = 20;
	else
		dim = 36;

	resultsFile << "Seq1\tSeq2\tChi-square\tdf\tp-value\tDelta_s\tDelta_ms\tSites\tStart\tEnd" << endl;
	for (unsigned int windowStart = 0; windowStart + windowSize <= cols; windowStart += windowStep)
	{
		vector<double> bowker_mat(n * n, 0);
		vector<double> ds_mat(n * n, 0);
		vector<double> dms_mat(n * n, 0);
		for (unsigned int k = 0; k < n; k++) // 1st sequence
		{
			for (unsigned int l = k + 1; l < n; l++) // 2nd sequence
			{
				unsigned int sum = 0;
				vector<int> dm(dim * dim, 0);
				for (unsigned int m = windowStart; m < windowStart + windowSize; m++)
				{
					Site *s = _sites[m];
					int c1 = s->getPos(k);
					int c2 = s->getPos(l);
					if (s->charIsUnambiguous(c1) && s->charIsUnambiguous(c2))
					{
						dm[c1 * dim + c2]++;
						sum++;
					}
				}

				unsigned int df = 0;
				double bowker = .0;
				for (int i = 0; i < dim; i++)
					for (int j = i + 1; j < dim; j++)
					{
						double dm_ij = dm[i * dim + j];
						double dm_ji = dm[j * dim + i];
						if ((dm[i * dim + j] + dm[j * dim + i]) > 0)
						{
							df++;
							bowker += ((dm_ij - dm_ji) * (dm_ij - dm_ji)) / (dm_ij + dm_ji);
						}
					}

				double Q = 1.0;
				if (df > 0)
					Q = gammq((df / 2.0), (bowker / 2.0));
				bowker_mat[k * n + l] = Q;
				bowker_mat[l * n + k] = Q;

				double delta_s = 0.0;
				for (int i = 0; i < dim; i++)
					for (int j = i + 1; j < dim; j++)
					{
						double x = dm[j * dim + i] - dm[i * dim + j];
						delta_s += (x / sum) * (x / sum);
					}
				delta_s = sqrt(delta_s);
				ds_mat[k * n + l] = delta_s;
				ds_mat[l * n + k] = delta_s;

				double delta_ms = 0.0;
				for (int i = 0; i < dim; i++)
				{
					double row = 0.0;
					double col = 0.0;
					for (int j = 0; j < dim; j++)
					{
						row += dm[i * dim + j];
						col += dm[j * dim + i];
					}
					delta_ms += ((row - col) / sum) * ((row - col) / sum);
				}
				delta_ms = sqrt(delta_ms) / sqrt(2.0);
				dms_mat[k * n + l] = delta_ms;
				dms_mat[l * n + k] = delta_ms;

				resultsFile << _alignment[k].getName() << "\t" << _alignment[l].getName() << "\t" << scientific << bowker << "\t" << df << "\t" << Q << "\t" << delta_s << "\t" << delta_ms << "\t" << sum << "\t" << windowStart << "\t" << windowStart
						+ windowSize - 1 << endl;
			}

		}

		if (extended)
		{
			bowkerFile.flags(ios::left);
			delta_sFile.flags(ios::left);
			delta_msFile.flags(ios::left);
			bowkerFile << windowStart << "-" << windowStart + windowSize - 1;
			delta_sFile << windowStart << "-" << windowStart + windowSize - 1;
			delta_msFile << windowStart << "-" << windowStart + windowSize - 1;

			for (unsigned int l = 0; l < n; l++)
			{
				bowkerFile << "\t" << setw(12) << _alignment[l].getName();
				delta_sFile << "\t" << setw(12) << _alignment[l].getName();
				delta_msFile << "\t" << setw(12) << _alignment[l].getName();
			}
			bowkerFile << endl;
			delta_sFile << endl;
			delta_msFile << endl;

			for (unsigned int k = 0; k < n; k++)
			{
				bowkerFile << setw(12) << _alignment[k].getName();
				delta_sFile << setw(12) << _alignment[k].getName();
				delta_msFile << setw(12) << _alignment[k].getName();
				for (unsigned int l = 0; l < n; l++)
				{
					bowkerFile << "\t" << scientific << bowker_mat[k * n + l];
					delta_sFile << "\t" << scientific << ds_mat[k * n + l];
					delta_msFile << "\t" << scientific << dms_mat[k * n + l];
				}
				bowkerFile << endl;
				delta_sFile << endl;
				delta_msFile << endl;
			}
		}
	}
	resultsFile.close();
	if (extended)
	{
		bowkerFile.close();
		delta_sFile.close();
		delta_msFile.close();
	}
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


void Alignment::computeBasicScores()
{
	cout << endl;
	cout << "Computing basic scores..." << endl;

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


void Alignment::computeCompatibilityScores(int randomizations)
{
	cout << endl;
	cout << "Computing compatibility scores, doing " << randomizations << " randomizations..." << endl;
	cout << "0%" << flush;

	long t1 = time(NULL);
	long lastTime = t1;
	unsigned long n = _informativeSites.size();
	unsigned long total = n * (n - 1) / 2 + n * n * randomizations;
	unsigned long count = 0;

	srand( t1 );
#ifdef _OPENMP
#pragma omp parallel for shared(count)
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
		count += n - i - 1;
		if (omp_get_thread_num() == 0)
		{
			long t2 = time(NULL);
			if (t2 > lastTime)
			{
				long elapsed = t2 - t1;
				long eta = (elapsed * total) / count - elapsed;
				cout << "\r" << count * 100 / total << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;

			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for shared(count)
#endif
	for (unsigned int i = 0; i < n; i++)
	{
		_informativeSites[i]->computeCo(n);

		if (randomizations)
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
				if (_informativeSites[i]->getComp() <= comp)
					poc++;
			}

			count += randomizations * n;
			if (omp_get_thread_num() == 0)
			{
				long t2 = time(NULL);
				if (t2 > lastTime)
				{
					long elapsed = t2 - t1;
					long eta = (elapsed * total) / count - elapsed;
					cout << "\r" << count * 100 / total << "%\tTime elapsed: " << printTime(elapsed) << "\tETA: " << printTime(eta) << "  " << flush;
				}
			}
			_informativeSites[i]->computePOC(poc, randomizations);
		}
	}

	long t2 = time(NULL);
	cout << "\rDone, taking " << printTime(t2-t1) << "                         " << endl;
}

Alignment Alignment::getModifiedAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy)
{
	Alignment a(_dataType);
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

	return a;
}

void Alignment::writeSummary(string prefix)
{
	string fileName = prefix + ".sites.csv";
	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("\n\nError, cannot open file " + fileName);
	cout << "Writing site summary to " << fileName << endl;

	set<int> bases;
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		BaseOccurenceMap f = s->getFrequencies();
		for (BaseOccurenceMapIterator it = f.begin(); it != f.end(); it++)
			bases.insert(it->first);
	}

	file << "Site No.,Smin,Entropy,OV,Co,Poc";
	Site* s = _informativeSites[0];
	for (set<int>::iterator it = bases.begin(); it != bases.end(); it++)
		file << ",f(" << s->mapNumToChar(*it) << ")";
	file << ",f(?)" << endl;

	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		BaseOccurenceMap f = s->getFrequencies();
		file << s->getCols()[0] + 1 << "," << s->getSmin() << "," << fixed << s->getEntropy() << "," << s->getOV() << "," << s->getCo() << "," << s->getPOC();
		for (set<int>::iterator it = bases.begin(); it != bases.end(); it++)
			file << "," << (((double) f[*it]) / getNumOfRows());
		file << "," << ((double) s->getAmbiguousCount()) / getNumOfRows() << endl;
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

	switch (type) {
		case 0:
			cout << "Writing Phylip alignment to " << fileName << endl;
			file << getNumOfRows() << " " << getNumOfCols() << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
				file << setw(10) << left << _alignment[i].getName() << _alignment[i].getSequence() << endl;
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
