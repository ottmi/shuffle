#ifndef AASITE_H_
#define AASITE_H_
#include "Site.h"

class AASite: public Site
{
public:
	AASite(vector<Sequence>* alignment, int col);
	AASite(vector<char> site, int col);
	virtual ~AASite();

private:
	char mapNumToChar(char c);
	char mapCharToNum(char c);
};

#endif /* AASITE_H_ */
