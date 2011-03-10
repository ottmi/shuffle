#ifndef AASITE_H_
#define AASITE_H_
#include "Site.h"
#include "Alignment.h"

class AASite: public Site
{
public:
	AASite(Alignment *alignment, int col);
	AASite(string site, int col);
	virtual ~AASite();
};

#endif /* AASITE_H_ */
