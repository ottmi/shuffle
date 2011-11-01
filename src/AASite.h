#ifndef AASITE_H_
#define AASITE_H_
#include "Site.h"

class AASite: public Site
{
public:
	AASite(vector<Sequence>* alignment, int offset, Options *options);
	AASite(vector<unsigned int> site);
	virtual ~AASite();

private:
	string mapNumToChar(unsigned int n);
	unsigned int mapCharToNum(string s);
};

#endif /* AASITE_H_ */
