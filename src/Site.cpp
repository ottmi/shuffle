#include <iostream>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "globals.h"
#include "AASite.h"
#include "DNASite.h"
#include "Sequence.h"
#include "Site.h"

Site::Site()
{
}


Site::~Site()
{
}


double factorial(int n)
{
	double m = (double) n;
	for (int i = n - 1; i > 1; i--)
		m *= i;

	return m;
}


void Site::initialize(vector<Sequence>* alignment)
{
	srand ( time(NULL) );

	_compSites = 0;
	_coScore = .0;

	_unambiguousCount = 0;
	for (unsigned int i = 0; i < alignment->size(); i++)
	{
		Sequence sequence = alignment->at(i);
		char c = this->mapCharToNum(sequence.getSequence().at(_col));

		if (charIsUnambiguous(c))
		{
			_r[c]++;
			_unambiguousCount++;
		}
		_site.push_back(c);
	}

	int informative = 0;
	for (char c=0; c<=_unambiguousThreshold; c++)
	{
		BaseOccurenceMapIterator it = _r.find(c);
		if (it != _r.end() && it->second >= 2)
			informative++;
	}

	if (informative >= 2)
		_isInformative = true;
	else
		_isInformative = false;
}


bool Site::checkCompatibility(Site* site)
{
	vector<char> s1 = _site;
	vector<char> s2 = site->getSite();

	list<pair<char, char> > pairs;
	list<pair<char, char> >::iterator it1;
	list<pair<char, char> >::iterator it2;

	for (unsigned int i = 0; i < s1.size(); i++)
	{
		if (charIsUnambiguous(s1[i]) && charIsUnambiguous(s2[i])) // contains no unknown, missing, or ambiguous characters
		{
			pair<char, char> p(s1[i], s2[i]);
			pairs.push_back(p);
		}
	}
/*
	pairs.sort();
	pairs.unique();
*/
	for (it1 = pairs.begin(); it1 != pairs.end(); it1++) // remove duplicates
	{
		it2 = it1;
		it2++;
		while (it2 != pairs.end())
		{
			if (it1->first == it2->first && it1->second == it2->second)
				it2 = pairs.erase(it2);
			else
				it2++;
		}
	}
/*
	for (it1 = pairs.begin(); it1 != pairs.end(); it1++)
		cerr << it1->first;
	cerr << endl;
	for (it1 = pairs.begin(); it1 != pairs.end(); it1++)
		cerr << it1->second;
	cerr << endl;
*/
	BaseOccurenceMap occ1;
	BaseOccurenceMap occ2;
	for (it1 = pairs.begin(); it1 != pairs.end(); it1++)
	{
		occ1[it1->first]++;
		occ2[it1->second]++;
	}

	bool changed = true;
	while (changed)
	{
		changed = false;
		for (it1 = pairs.begin(); it1 != pairs.end(); it1++)
		{
			if (occ1[it1->first] == 1)
			{
				occ1[it1->first]--;
				occ2[it1->second]--;
				it1 = pairs.erase(it1);
				changed = true;
			}
		}

		for (it2 = pairs.begin(); it2 != pairs.end(); it2++)
		{
			if (occ2[it2->second] == 1)
			{
				occ1[it2->first]--;
				occ2[it2->second]--;
				it2 = pairs.erase(it2);
				changed = true;
			}
		}
	}

	return (pairs.size() == 0);
}


void Site::incComp()
{
#ifdef _OPENMP
	#pragma omp atomic
#endif
	_compSites++;
}


void Site::setPOC(double poc)
{
	_poc = poc;
	if (verbose)
		cout << "Site #" << _col << ": compSites=" << _compSites << " co=" << _coScore << " poc=" << _poc << " entropy=" << _entropy << " smin=" << _smin << endl;
}


void Site::computeScores(unsigned int cols)
{
	_smin = _r.size() - 1;

	BaseOccurenceMapIterator r_it;
	map<int, int> t;
	map<int, int>::const_iterator t_it;

	double prod_r = 1.0;
	for (r_it = _r.begin(); r_it != _r.end(); r_it++)
	{
		t[r_it->second]++;
		prod_r *= factorial(r_it->second);
	}

	double prod_t = 1.0;
	for (t_it = t.begin(); t_it != t.end(); t_it++)
	{
		prod_t *= factorial(t_it->second);
	}
	_entropy = log(factorial(_unambiguousCount) / (prod_r * prod_t));

	unsigned int n = 0;
	for (char c=0; c<=_unambiguousThreshold; c++)
		n+= _r[c];

	unsigned int d = 0;
	for (char c=0; c<=_unambiguousThreshold; c++)
	{
		d+= (n - _r[c]) * _r[c];
	}
	double k = n*(n-1);
	_ov = (double) d/k;

	_coScore = ((double) _compSites) / cols;
}


Site* Site::randomize()
{
	vector<char> r;
	Site* randomizedSite = NULL;
	vector<char> positions = _site;

	for (unsigned int i=0; i<_site.size(); i++)
	{
		unsigned int j = rand() % positions.size();
		r.push_back(positions[j]);
		positions.erase(positions.begin()+j);
	}

	switch (_type)
	{
		case _DNA_DATA:
			randomizedSite = new DNASite(r, -1);
			break;
		case _AA_DATA:
			randomizedSite = new AASite(r, -1);
			break;
		default:
			cerr << "Unknown data type " << _type << " at Site::randomize()" << endl;
			break;
	}

	return randomizedSite;
}


bool Site::charIsUnambiguous(char c)
{
	return (c <= _unambiguousThreshold);
}
