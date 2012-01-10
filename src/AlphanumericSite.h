#ifndef ALPHANUMERICSITE_H_
#define ALPHANUMERICSITE_H_

#include "Site.h"

class AlphanumericSite: public Site
{
public:
	AlphanumericSite(int offset, Options *options);
	AlphanumericSite(vector<unsigned int> site);
	virtual ~AlphanumericSite();

private:
	string mapNumToChar(unsigned int n);
	unsigned int mapCharToNum(string s);
};

#endif /* ALPHANUMERICSITE_H_ */
