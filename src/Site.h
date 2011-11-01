#ifndef SITE_H_
#define SITE_H_

#include <string>
#include <map>
#include <vector>
#include "globals.h"
#include "Sequence.h"

using namespace std;

typedef map<unsigned int,int> BaseOccurenceMap;
typedef map<unsigned int,int>::iterator BaseOccurenceMapIterator;

class Site
{
public:
	Site();
	virtual ~Site();
	void initialize(vector<Sequence>* alignment, Options *options);
	void remove(unsigned int i);
	bool checkInformative();
	bool checkCompatibility(Site* site);
	void incComp();
	void computeScores(unsigned int cols);
	void computeCo(unsigned int cols);
	void computePOC(int poc, int randomizations);
	bool compare(Site* s);
	Site* randomize();
	void setPOC(double poc);
	vector<unsigned int> getSite() { return _site; };
	vector<int> getCols() { return _cols; };
	unsigned int getPos(unsigned int pos) { return _site[pos]; };
	bool isInformative() { return _isInformative; };
	int getComp() { return _compSites; };
	double getCo() { return _coScore; };
	double getPOC() { return _poc; };
	double getEntropy() { return _entropy; };
	int getSmin() { return _smin; };
	double getOV() { return _ov; };
	int getUnambiguousCount() { return _unambiguousCount; };
	int getAmbiguousCount() { return _ambiguousCount; };
	BaseOccurenceMap& getFrequencies() { return _r; };
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
	BaseOccurenceMap _r;
	int _unambiguousCount;
	int _ambiguousCount;
	vector<int> _cols;
	int _compSites;
	bool _isInformative;
	double _coScore;
	double _poc;
	double _entropy;
	int _smin;
	double _ov;
};

#endif /* SITE_H_ */
