#ifndef DNASITE_H_
#define DNASITE_H_
#include "Site.h"

class DNASite: public Site
{
public:
	DNASite(vector<Sequence>* alignment, vector<int> grouping, int offset);
	DNASite(vector<int> site);
	virtual ~DNASite();

private:
	string mapNumToChar(int n);
	int mapCharToNum(string s);
};

#endif /* DNASITE_H_ */
