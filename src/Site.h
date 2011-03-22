#ifndef SITE_H_
#define SITE_H_

#include <string>
#include <map>
#include <vector>
#include "Sequence.h"

using namespace std;

typedef map<char,int> BaseOccurenceMap;
typedef map<char,int>::iterator BaseOccurenceMapIterator;

class Site
{
public:
	Site();
	virtual ~Site();
	void initialize(vector<Sequence>* alignment);
	BaseOccurenceMap getBaseOccurences();
	bool checkCompatibility(Site* site);
	void incComp();
	void computeCompScore(unsigned int cols);
	Site* randomize();
	void setPOC(double poc);
	virtual bool charIsUnambiguous(char c);
	string getString() { return _site; };
	int getCol() { return _col; };
	bool isInformative() { return _isInformative; };
	int getComp() { return _compSites; };
	double getCo() { return _coScore; };
	double getPOC() { return _poc; };
	double getEntropy() { return _entropy; };
	int getSmin() { return _smin; };
	double getOV() { return _ov; };

protected:
	int _type; // 0=DNA, 1=AA
	string _unambiguousCharacters;
	string _ambiguousCharacters;
	string _missingCharacters;
	int _col;
	string _site;
	int _compSites;
	bool _isInformative;
	double _coScore;
	double _poc;
	double _entropy;
	int _smin;
	double _ov;
};

#endif /* SITE_H_ */
