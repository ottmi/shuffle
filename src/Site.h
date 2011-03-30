#ifndef SITE_H_
#define SITE_H_

#include <string>
#include <map>
#include <vector>
#include "Sequence.h"

using namespace std;

typedef map<int,int> BaseOccurenceMap;
typedef map<int,int>::iterator BaseOccurenceMapIterator;

class Site
{
public:
	Site();
	virtual ~Site();
	void initialize(vector<Sequence>* alignment);
	BaseOccurenceMap getBaseOccurences();
	bool checkCompatibility(Site* site);
	void incComp();
	void computeScores(unsigned int cols);
	Site* randomize();
	void setPOC(double poc);
	vector<int> getSite() { return _site; };
	vector<int> getCols() { return _cols; };
	bool isInformative() { return _isInformative; };
	int getComp() { return _compSites; };
	double getCo() { return _coScore; };
	double getPOC() { return _poc; };
	double getEntropy() { return _entropy; };
	int getSmin() { return _smin; };
	double getOV() { return _ov; };
	bool charIsUnambiguous(int n);
	virtual string mapNumToChar(int n) =0;
	virtual int mapCharToNum(string s) =0;
	string toString();
	string toNumString();
	string colsToString();

protected:
	int _type; // 0=DNA, 1=AA
	vector<int> _site;
	char _unambiguousThreshold;
	BaseOccurenceMap _r;
	int _unambiguousCount;
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
