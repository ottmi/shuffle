#include <iostream>
#include <iomanip>
#include <set>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include "globals.h"
#include "helper.h"
#include "AASite.h"
#include "DNASite.h"
#include "AlphanumericSite.h"
#include "Sequence.h"
#include "Site.h"
#ifdef _MPFR
#include <mpfr.h>
#endif


Site::Site()
{
}


Site::~Site()
{
}

/* This function is called from Alignment::recv() and Site::randomize() via
 * AASite::AASite(), AlphanumericSite::AlphanumericSite(), and DNASite::DNASite() */
void Site::initialize()
{
#ifdef _MPI
	_compatibleSites = 0;
#endif
	_unambiguousCount = 0;
	_ambiguousCount = 0;
	_coScore = .0;
	_poc = .0;
	_entropy = .0;
	_smin = 0;
	_ov = .0;
	createRandomSeed(_randomSeed, _cols[0]);
}

/* This function is called from Alignment::collectSites() */
void Site::initialize(vector<Sequence>* alignment)
{
#ifdef _MPI
	_compatibleSites = 0;
#endif
	_unambiguousCount = 0;
	_ambiguousCount = 0;
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
	createRandomSeed(_randomSeed, _cols[0]);
}


void Site::remove(unsigned int i)
{
	_site.erase(_site.begin() + i);
}


bool Site::checkInformative()
{
	_unambiguousCount = 0;
	_baseOccurences.clear();
	for (unsigned int i = 0; i < _site.size(); i++)
	{
		unsigned int c = _site[i];
		_pattern[c].insert(i);

		if (charIsUnambiguous(c))
		{
			_baseOccurences[c]++;
			_unambiguousCount++;
		}
	}
	_ambiguousCount = _site.size()-_unambiguousCount;

	int informative = 0;
	for (BaseOccurenceMapIterator it=_baseOccurences.begin(); it!=_baseOccurences.end(); it++)
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

	if (verbose >= 5)
		cout << "   " << toString() << endl;
	if (verbose >= 6)
		cout << "   " << toNumString() << endl;

	return _isInformative;
}


unsigned int getFirst(unsigned long long x)
{
	return x >> 32;
}

unsigned int getSecond(unsigned long long x)
{
	return x & 0xffffffff;
}

bool Site::checkCompatibility(Site* site)
{
	vector<unsigned int> s1 = _site;
	vector<unsigned int> s2 = site->getSite();

	set<unsigned long long> pairs;
	for (unsigned int i = 0; i < s1.size(); i++)
	{
		if (charIsUnambiguous(s1[i]) && charIsUnambiguous(s2[i])) // contains no unknown, missing, or ambiguous characters
		{
			unsigned long long p = s1[i];
			p = (p << 32) | s2[i];
			pairs.insert(p);
		}
	}

	BaseOccurenceMap occ1;
	BaseOccurenceMap occ2;
	for (set<unsigned long long>::iterator it = pairs.begin(); it != pairs.end(); it++)  // count occurences
	{
		occ1[getFirst(*it)]++;
		occ2[getSecond(*it)]++;
	}

	bool changed = true;
	while (changed) // remove per-site unique characters
	{
		changed = false;
		set<unsigned long long>::iterator it = pairs.begin();
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


#ifdef _MPI
void Site::addCompatibleSite(int site)
{
	_compatibleSites++;
}

void Site::removeCompatibleSite(int site)
{
	_compatibleSites--;
}

int Site::getComp()
{
    return _compatibleSites;
}

void Site::setComp(unsigned int val)
{
    _compatibleSites = val;
}
#else
void Site::addCompatibleSite(int site)
{
#ifdef _OPENMP
#pragma omp critical
#endif
	_compatibleSites.insert(site);
}

void Site::removeCompatibleSite(int site)
{
	_compatibleSites.erase(site);
}

int Site::getComp()
{
    return _compatibleSites.size();

}
#endif


double Site::checkPattern(Site* site)
{
    /*
     * Carla A. Cummins and James O. McInerney:
     * "A Method for Inferring the Rate of Evolution of Homologous Characters that Can Potentially
     * Improve Phylogenetic Inference, Resolve Deep Divergence and Correct Systematic Biases"
     * Syst. Biol. 60(6):833–844, 2011
     */

    SitePattern p_i = getPattern();
    SitePattern p_j = site->getPattern();
    set<unsigned int>::iterator base;
    unsigned int count = 0;

    if (verbose >= 5)
    {
	cout << "  Site [" << getCols()[0] << "] ";
	for (SitePatternIterator part = p_i.begin(); part != p_i.end(); part++)
	{
	    for (base = part->second.begin(); base != part->second.end(); base++)
		cout << *base << " ";
	    cout << "   ";
	}
	cout << endl;
	cout << "  Site [" << site->getCols()[0] << "] ";
	for (SitePatternIterator part = p_j.begin(); part != p_j.end(); part++)
	{
	    for (base = part->second.begin(); base != part->second.end(); base++)
		cout << *base << " ";
	    cout << "   ";
	}
	cout << endl;
    }

    for (SitePatternIterator part1 = p_j.begin(); part1 != p_j.end(); part1++)
    {
	for (SitePatternIterator part2 = p_i.begin(); part2 != p_i.end(); part2++)
	{
	    unsigned int matches = 0;
	    for (base = part1->second.begin(); base != part1->second.end(); base++)
	    {
		if (part2->second.count(*base))
		    matches++;
		else
		    break;
	    }

	    if (matches == part1->second.size())
		count++;
	}
    }

    double pa_ij = ((double) count) / p_j.size();
    if (verbose >= 5)
	cout << "  patterns=" << p_j.size() << " count=" << count << " pa_ij=" << pa_ij << endl;

    return pa_ij;
}


void Site::computeScores(unsigned int cols)
{
	_smin = _baseOccurences.size() - 1;
	map<unsigned int, int> t;
	map<unsigned int, int>::const_iterator t_it;

#ifdef _MPFR
	mpfr_set_default_prec(512);
	mpfr_t tmp1, tmp2;
	mpfr_init(tmp1);
	mpfr_init(tmp2);
	mpfr_t prod_r;
	mpfr_init_set_ui(prod_r, 1, GMP_RNDN);			/* prod_r = 1 */
	for (BaseOccurenceMapIterator it = _baseOccurences.begin(); it != _baseOccurences.end(); it++)
	{
		t[it->second]++;
		mpfr_fac_ui(tmp1, it->second, GMP_RNDN);	/* tmp1 = factorial(it->second) */
		mpfr_set(tmp2, prod_r, GMP_RNDN);		/* tmp2 = prod_r */
		mpfr_mul(prod_r, tmp1, tmp2, GMP_RNDN);		/* prod_r = tmp1 * tmp2 */
	}

	mpfr_t prod_t;
	mpfr_init_set_ui(prod_t, 1, GMP_RNDN);			/* prod_t = 1 */
	for (t_it = t.begin(); t_it != t.end(); t_it++)
	{
		mpfr_fac_ui(tmp1, t_it->second, GMP_RNDN);	/* tmp1 = factorial(t_it->second) */
		mpfr_set(tmp2, prod_t, GMP_RNDN);		/* tmp2 = prod_t */
		mpfr_mul(prod_t, tmp1, tmp2, GMP_RNDN);		/* prod_t = tmp1 * tmp2 */
	}

	mpfr_t unamb;
	mpfr_init(unamb);
	mpfr_fac_ui(unamb, _unambiguousCount, GMP_RNDN);	/* unamb = factorial(_unambiguousCount) */

	mpfr_div(tmp1, unamb, prod_r, GMP_RNDN);		/* tmp1 = unamb / prod_r */
	mpfr_div(tmp2, tmp1, prod_t, GMP_RNDN);			/* tmp2 = tmp1 / prod_t */
	mpfr_log(tmp1, tmp2, GMP_RNDN);				/* tmp1 = log(tmp2) */

	_entropy = mpfr_get_d(tmp1, GMP_RNDN);
	mpfr_clears(tmp1, tmp2, prod_r, prod_t, unamb, (mpfr_ptr) 0);
#else
	double prod_r = 1;
	for (BaseOccurenceMapIterator it = _baseOccurences.begin(); it != _baseOccurences.end(); it++)
	{
		t[it->second]++;
		prod_r *= factorial(it->second);
	}

	double prod_t = 1.0;
	for (t_it = t.begin(); t_it != t.end(); t_it++)
	{
		prod_t *= factorial(t_it->second);
	}
	double unamb = factorial(_unambiguousCount);
	double div =  (unamb / prod_r) / prod_t;
	_entropy = log(div);
#endif

	unsigned int n = 0;
	for (BaseOccurenceMapIterator it = _baseOccurences.begin(); it != _baseOccurences.end(); it++)
		n+= it->second;

	unsigned int d = 0;
	for (BaseOccurenceMapIterator it=_baseOccurences.begin(); it!=_baseOccurences.end(); it++)
		if (charIsUnambiguous(it->first))
			d+= (n - it->second) * it->second;

	double k = n*(n-1);
	_ov = (double) d/k;
}


void Site::computeCo(unsigned int cols)
{
	_coScore = ((double) getComp()) / (cols-1);
	if (verbose >= 5)
		cout << "\rSite [" << colsToString() << "]: compSites=" << getComp() << " co=" << _coScore << endl;
}


void Site::computePOC(int poc, int randomizations)
{
	_poc = ((double) poc) / randomizations;
	if (verbose >= 5)
		cout << "\rSite [" << colsToString() << "]: poc=" << _poc << endl;
}


Site* Site::randomize()
{
	vector<unsigned int> r(_site.size(), 0);
	Site* randomizedSite = NULL;
	vector<unsigned int> positions = _site;

	for (unsigned int i=0; i<_site.size(); i++)
	{
		unsigned int j = nrand48(_randomSeed) % positions.size();
		r[i] = positions[j];
		positions.erase(positions.begin()+j);
	}

	switch (_type)
	{
		case _DNA_DATA:
			randomizedSite = new DNASite(-1, r);
			break;
		case _AA_DATA:
			randomizedSite = new AASite(-1, r);
			break;
		case _ALPHANUM_DATA:
			randomizedSite = new AlphanumericSite(-1, r);
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
		if ((char) (n & 255) > _unambiguousThreshold)
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
