/*
 * DNASite.h
 *
 *  Created on: Mar 8, 2011
 *      Author: ottmi
 */

#ifndef DNASITE_H_
#define DNASITE_H_
#include "Alignment.h"
#include "Site.h"

class DNASite: public Site
{
public:
	DNASite(Alignment *alignment, int col);
	DNASite(string site, int col);
	virtual ~DNASite();
};

#endif /* DNASITE_H_ */
