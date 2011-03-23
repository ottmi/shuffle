#ifndef DNASITE_H_
#define DNASITE_H_
#include "Site.h"

class DNASite: public Site
{
public:
	DNASite(vector<Sequence>* alignment, int col);
	DNASite(vector<char> site, int col);
	virtual ~DNASite();

private:
	char mapNumToChar(char c);
	char mapCharToNum(char c);
};

#endif /* DNASITE_H_ */
