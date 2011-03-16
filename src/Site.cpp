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


long factorial(long n)
{
	for (long i = n - 1; i > 1; i--)
		n *= i;
	return n;
}


void Site::initialize(vector<Sequence>* alignment)
{
	srand ( time(NULL) );

	_compSites = 0;
	_compScore = .0;

	BaseOccurenceMap r;
	BaseOccurenceMapIterator r_it;
	map<int, int> t;
	map<int, int>::const_iterator t_it;

	for (unsigned int i = 0; i < alignment->size(); i++)
	{
		Sequence sequence = alignment->at(i);
		char c = sequence.getSequence().at(_col);
		r[c]++;
		_site += c;
	}
	_mnic = r.size() - 1;

	long prod_r = 1;
	for (r_it = r.begin(); r_it != r.end(); r_it++)
	{
		t[r_it->second]++;
		prod_r *= factorial(r_it->second);
	}

	long prod_t = 1;
	for (t_it = t.begin(); t_it != t.end(); t_it++)
	{
		prod_t *= factorial(t_it->second);
	}
	_entropy = log((double) factorial(alignment->size()) / (prod_r * prod_t));

	int informative = 0;
	for (unsigned int i=0; i<_unambiguousCharacters.length(); i++)
	{
		char c = _unambiguousCharacters[i];
		if (r[c] >= 2)
			informative++;
	}
	if (informative >= 2)
		_isInformative = true;
	else
		_isInformative = false;
}


bool Site::checkCompatibility(Site* site)
{
	string s1 = _site;
	string s2 = site->getString();

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
	_compSites++;
}


void Site::setPOC(double poc)
{
	_poc = poc;
	if (verbose)
		cout << "Site #" << _col << ": compSites=" << _compSites << " co=" << _compScore << " poc=" << _poc << " entropy=" << _entropy << " mnic=" << _mnic << endl;
}


void Site::computeCompScore(unsigned int cols)
{
	_compScore = ((double) _compSites) / cols;
}


Site* Site::randomize()
{
	string r;
	Site* randomizedSite = NULL;

	vector <char> positions;

	for (unsigned int i=0; i<_site.length(); i++)
		positions.push_back(_site[i]);

	for (unsigned int i=0; i<_site.length(); i++)
	{
		unsigned int j = rand() % positions.size();
		r+= positions[j];
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
	return (_unambiguousCharacters.find(c) != string::npos);
}
