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
	double getCompScore() { return _compScore; };
	double getPOC() { return _poc; };
	double getEntropy() { return _entropy; };
	int getMNIC() { return _mnic; };

protected:
	int _type; // 0=DNA, 1=AA
	string _unambiguousCharacters;
	string _ambiguousCharacters;
	string _missingCharacters;
	int _col;
	string _site;
	int _compSites;
	bool _isInformative;
	double _compScore;
	double _poc;
	double _entropy;
	int _mnic;
};

#endif /* SITE_H_ */
