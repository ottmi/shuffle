#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "globals.h"
#include "helper.h"
#include "AASite.h"
#include "DNASite.h"
#include "AlphanumericSite.h"
#include "Sequence.h"
#include "Site.h"

Site::Site()
{
}


Site::~Site()
{
}


void Site::initialize()
{
	_unambiguousCount = 0;
	_ambiguousCount = 0;
	_compSites = 0;
	_coScore = .0;
	_poc = .0;
	_entropy = .0;
	_smin = 0;
	_ov = .0;
}


void Site::initialize(vector<Sequence>* alignment)
{
	_unambiguousCount = 0;
	_ambiguousCount = 0;
	_compSites = 0;
	_coScore = .0;
	_poc = .0;
	_entropy = .0;
	_smin = 0;
	_ov = .0;

	for (unsigned int i = 0; i < alignment->size(); i++)
	{
		Sequence sequence = alignment->at(i);
		int c = this->mapCharToNum(sequence.getColumns(_cols));
		_site.push_back(c);
	}
}

void Site::remove(unsigned int i)
{
	_site.erase(_site.begin() + i);
}


bool Site::checkInformative()
{
	_unambiguousCount = 0;
	_r.clear();
	for (unsigned int i = 0; i < _site.size(); i++)
	{
		int c = _site[i];

		if (charIsUnambiguous(c))
		{
			_r[c]++;
			_unambiguousCount++;
		}
	}
	_ambiguousCount = _site.size()-_unambiguousCount;

	int informative = 0;
	for (BaseOccurenceMapIterator it=_r.begin(); it!=_r.end(); it++)
		if (charIsUnambiguous(it->first) && it->second >= 2)
			informative++;

	if (informative >= 2)
	{
		if (verbose >= 2)
			cout << "Site [" << colsToString() << "] is informative" << endl;
		_isInformative = true;
	}
	else
	{
		if (verbose >= 3)
			cout << "Site [" << colsToString() << "] is uninformative" << endl;
		_isInformative = false;
	}

	if (verbose >= 3)
		cout << "   " << toString() << endl;
	if (verbose >= 4)
		cout << "   " << toNumString() << endl;

	return _isInformative;
}


int getFirst(unsigned long x)
{
	return x >> 32;
}

int getSecond(unsigned long x)
{
	return x & 0xffffffff;
}

bool Site::checkCompatibility(Site* site)
{
	vector<int> s1 = _site;
	vector<int> s2 = site->getSite();

	set<unsigned long> pairs;
	set<unsigned long>::iterator it1;
	set<unsigned long>::iterator it2;

	for (unsigned int i = 0; i < s1.size(); i++)
	{
		if (charIsUnambiguous(s1[i]) && charIsUnambiguous(s2[i])) // contains no unknown, missing, or ambiguous characters
		{
			unsigned long p = s1[i];
			p = (p << 32) | s2[i];
			pairs.insert(p);
		}
	}

	BaseOccurenceMap occ1;
	BaseOccurenceMap occ2;
	for (it1 = pairs.begin(); it1 != pairs.end(); it1++)  // count occurences
	{
		occ1[getFirst(*it1)]++;
		occ2[getSecond(*it1)]++;
	}

	bool changed = true;
	while (changed) // remove per-site unique characters
	{
		changed = false;
		for (it1 = pairs.begin(); it1 != pairs.end(); it1++)
		{
			int first = getFirst(*it1);
			if (occ1[first] == 1)
			{
				int second = getSecond(*it1);
				occ1[first]--;
				occ2[second]--;
				pairs.erase(it1);
				changed = true;
			}
		}

		for (it2 = pairs.begin(); it2 != pairs.end(); it2++)
		{
			int second = getSecond(*it2);
			if (occ2[second] == 1)
			{
				int first = getFirst(*it2);
				occ1[first]--;
				occ2[second]--;
				pairs.erase(it2);
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


void Site::computeScores(unsigned int cols)
{
	_smin = _r.size() - 1;
	map<int, int> t;
	map<int, int>::const_iterator t_it;

	unsigned int n = 0;
	bignum prod_r = 1;
	for (BaseOccurenceMapIterator it = _r.begin(); it != _r.end(); it++)
	{
		n+= it->second;
		t[it->second]++;
		prod_r *= factorial(it->second);
	}

	bignum prod_t = 1.0;
	for (t_it = t.begin(); t_it != t.end(); t_it++)
	{
		prod_t *= factorial(t_it->second);
	}
	bignum unamb = factorial(_unambiguousCount);
	bignum div =  (unamb / prod_r) / prod_t;
	_entropy = bignum_log(div);

	unsigned int d = 0;
	for (BaseOccurenceMapIterator it=_r.begin(); it!=_r.end(); it++)
		if (charIsUnambiguous(it->first))
			d+= (n - it->second) * it->second;

	double k = n*(n-1);
	_ov = (double) d/k;
}


void Site::computeCo(unsigned int cols)
{
	_coScore = ((double) _compSites) / (cols-1);
}


void Site::computePOC(int poc, int randomizations)
{
	_poc = ((double) poc) / randomizations;
	if (verbose >= 5)
		cout << "\rSite [" << colsToString() << "]: compSites=" << _compSites << " co=" << _coScore << " poc=" << _poc << " entropy=" << _entropy << " smin=" << _smin << endl;
}


Site* Site::randomize()
{
	vector<int> r;
	Site* randomizedSite = NULL;
	vector<int> positions = _site;

	for (unsigned int i=0; i<_site.size(); i++)
	{
		unsigned int j = rand() % positions.size();
		r.push_back(positions[j]);
		positions.erase(positions.begin()+j);
	}

	switch (_type)
	{
		case _DNA_DATA:
			randomizedSite = new DNASite(r);
			break;
		case _AA_DATA:
			randomizedSite = new AASite(r);
			break;
		case _ALPHANUM_DATA:
			randomizedSite = new AlphanumericSite(r);
			break;
		default:
			cerr << "Unknown data type " << _type << " at Site::randomize()" << endl;
			break;
	}

	return randomizedSite;
}


bool Site::compare(Site* s)
{
	vector<int> s1 = s->getSite();
	vector<int> s2 = getSite();

	if (s1.size() != s2.size())
		return false;

	for (unsigned int i = 0; i < s1.size(); i++)
	{
		if (s1[i] != s2[i])
			return false;
	}

	return true;
}


bool Site::charIsUnambiguous(int n)
{
	for (unsigned int i = 0; i < _cols.size(); i++)
	{
		if ((n & 255) > _unambiguousThreshold)
			return false;
		n = n >> 8;
	}
	return true;
}


string Site::toString()
{
	string s;

	for (unsigned int i = 0; i < _site.size()-1; i++)
	{
		s+= mapNumToChar(_site[i]);
		s+= ", ";
	}
	s+= mapNumToChar(_site[_site.size()-1]);

	return s;
}


string Site::toNumString()
{
	stringstream s;

	s <<  hex << setfill('0');
	for (unsigned j=0; j<_site.size()-1; j++)
		s << setw(_cols.size()*2) << _site[j] << ", ";
	s << setw(_cols.size()*2) << _site[_site.size()-1];

	return s.str();
}


string Site::colsToString()
{
	stringstream s;

	for (unsigned int i = 0; i < _cols.size()-1; i++)
		s << _cols[i] << ",";
	s << _cols[_cols.size()-1];

	return s.str();
}
