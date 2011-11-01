#ifndef DNASITE_H_
#define DNASITE_H_
#include "Site.h"

class DNASite: public Site
{
public:
	DNASite(vector<Sequence>* alignment, int offset, Options *options);
	DNASite(vector<unsigned int> site);
	virtual ~DNASite();

private:
	string mapNumToChar(unsigned int n);
	unsigned int mapCharToNum(string s);
};

#endif /* DNASITE_H_ */
