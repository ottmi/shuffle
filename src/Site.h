#ifndef SITE_H_
#define SITE_H_

#include <string>
#include <map>
#include <vector>
#include <set>
#include "globals.h"
#include "Sequence.h"

using namespace std;

typedef map<unsigned int,int> BaseOccurenceMap;
typedef map<unsigned int,int>::iterator BaseOccurenceMapIterator;

typedef map<unsigned int, set<unsigned int> > SitePattern;
typedef map<unsigned int, set<unsigned int> >::iterator SitePatternIterator;

class Site
{
public:
	Site();
	virtual ~Site();
	void initialize();
	void initialize(vector<Sequence>* alignment);
	void remove(unsigned int i);
	bool checkInformative();
	bool checkCompatibility(Site* site);
	void addCompatibleSite(int site);
	void removeCompatibleSite(int site);
	double checkPattern(Site* site);
	void computeScores(unsigned int cols);
	void computeCo(unsigned int cols);
	void computePOC(int poc, int randomizations);
	bool compare(Site* s);
	Site* randomize();
	void setPOC(double poc) { _poc = poc; };
	void setR(double r) { _r = r; };
	double getR() { return _r; };
	vector<unsigned int> getSite() { return _site; };
	vector<int> getCols() { return _cols; };
	unsigned int getPos(unsigned int pos) { return _site[pos]; };
	bool isInformative() { return _isInformative; };
	int getComp();
	void setComp(unsigned int val);
	double getCo() { return _coScore; };
	double getPOC() { return _poc; };
	double getEntropy() { return _entropy; };
	int getSmin() { return _smin; };
	double getOV() { return _ov; };
	void addRandomizedCo(double co) { _randomizedCo.push_back(co); };
	vector<double>& getRandomizedCo() { return _randomizedCo; };
	SitePattern& getPattern() { return _pattern; };
	int getUnambiguousCount() { return _unambiguousCount; };
	int getAmbiguousCount() { return _ambiguousCount; };
	BaseOccurenceMap& getFrequencies() { return _baseOccurences; };
	bool charIsUnambiguous(unsigned int n);
	virtual string mapNumToChar(unsigned int n) =0;
	virtual unsigned int mapCharToNum(string s) =0;
	string toString();
	string toNumString();
	string colsToString();

protected:
	int _type; // 0=DNA, 1=AA
	vector<unsigned int> _site;
	char _unambiguousThreshold;
	BaseOccurenceMap _baseOccurences;
	SitePattern _pattern;
	int _unambiguousCount;
	int _ambiguousCount;
	vector<int> _cols;
#ifdef _MPI
	unsigned int _compatibleSites;
#else
	set<int> _compatibleSites;
#endif
	bool _isInformative;
	double _coScore;
	double _poc;
	double _entropy;
	int _smin;
	double _ov;
	double _r;
	vector<double> _randomizedCo;
	unsigned short _randomSeed[3];

};

#endif /* SITE_H_ */
