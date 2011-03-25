#ifndef AASITE_H_
#define AASITE_H_
#include "Site.h"

class AASite: public Site
{
public:
	AASite(vector<Sequence>* alignment, vector<int> grouping, int offset);
	AASite(vector<int> site);
	virtual ~AASite();

private:
	string mapNumToChar(int n);
	int mapCharToNum(string s);
};

#endif /* AASITE_H_ */
