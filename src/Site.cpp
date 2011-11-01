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


void Site::initialize(vector<Sequence>* alignment, Options *options)
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
		unsigned int c = this->mapCharToNum(sequence.getColumns(_cols));
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
		unsigned int c = _site[i];

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


unsigned int getFirst(unsigned long x)
{
	return x >> 32;
}

unsigned int getSecond(unsigned long x)
{
	return x & 0xffffffff;
}

bool Site::checkCompatibility(Site* site)
{
	vector<unsigned int> s1 = _site;
	vector<unsigned int> s2 = site->getSite();

	set<unsigned long> pairs;
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
	for (set<unsigned long>::iterator it = pairs.begin(); it != pairs.end(); it++)  // count occurences
	{
		occ1[getFirst(*it)]++;
		occ2[getSecond(*it)]++;
	}

	bool changed = true;
	while (changed) // remove per-site unique characters
	{
		changed = false;
		set<unsigned long>::iterator it = pairs.begin();
		while (it != pairs.end())
		{
			unsigned int first = getFirst(*it);
			if (occ1[first] == 1)
			{
				unsigned int second = getSecond(*it);
				occ1[first]--;
				occ2[second]--;
				pairs.erase(it++);
				changed = true;
			} else
			{
				it++;
			}
		}

		it = pairs.begin();
		while (it != pairs.end())
		{
			unsigned int second = getSecond(*it);
			if (occ2[second] == 1)
			{
				unsigned int first = getFirst(*it);
				occ1[first]--;
				occ2[second]--;
				pairs.erase(it++);
				changed = true;
			} else
			{
				it++;
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
	map<unsigned int, int> t;
	map<unsigned int, int>::const_iterator t_it;

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
	vector<unsigned int> r(_site.size(), 0);
	Site* randomizedSite = NULL;
	vector<unsigned int> positions = _site;

	for (unsigned int i=0; i<_site.size(); i++)
	{
		unsigned int j = rand() % positions.size();
		r[i] = positions[j];
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
	vector<unsigned int> s1 = s->getSite();
	vector<unsigned int> s2 = getSite();

	if (s1.size() != s2.size())
		return false;

	for (unsigned int i = 0; i < s1.size(); i++)
	{
		if (s1[i] != s2[i])
			return false;
	}

	return true;
}


bool Site::charIsUnambiguous(unsigned int n)
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
