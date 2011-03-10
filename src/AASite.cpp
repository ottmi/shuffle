#include "AASite.h"

AASite::AASite(Alignment *alignment, int col)
{
	_unambiguousCharacters = "ACDEFGHIKLMNPQRSTVWY";
	_ambiguousCharacters = "BJXZ?";
	_missingCharacters = '-';

	_type = 1;
	_col = col;

	initialize(alignment->getAlignment());
}

AASite::AASite(string site, int col)
{
	_unambiguousCharacters = "ACDEFGHIKLMNPQRSTVWY";
	_ambiguousCharacters = "BJXZ?";
	_missingCharacters = '-';

	_type = 1;
	_site = site;
	_col = col;
}

AASite::~AASite()
{
}
