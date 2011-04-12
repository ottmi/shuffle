#ifndef ALPHANUMERICSITE_H_
#define ALPHANUMERICSITE_H_

#include "Site.h"

class AlphanumericSite: public Site
{
public:
	AlphanumericSite(vector<Sequence>* alignment, vector<int> grouping, int offset);
	AlphanumericSite(vector<int> site);
	virtual ~AlphanumericSite();

private:
	string mapNumToChar(int n);
	int mapCharToNum(string s);
};

#endif /* ALPHANUMERICSITE_H_ */
