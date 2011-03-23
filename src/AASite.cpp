#include "AASite.h"

AASite::AASite(vector<Sequence>* alignment, int col)
{
/*
	_unambiguousCharacters = "ACDEFGHIKLMNPQRSTVWY";
	_ambiguousCharacters = "BJXZ?";
	_missingCharacters = '-';
	_unambiguousThreshold = 19;
*/
	_unambiguousThreshold = 19;
	_type = 1;
	_col = col;
	initialize(alignment);
}


AASite::AASite(vector<char> site, int col)
{
	_unambiguousThreshold = 19;
	_type = 1;
	_col = col;
	_site = site;
}


AASite::~AASite()
{
}


char AASite::mapNumToChar(char c)
{
	string map = "ACDEFGHIKLMNPQRSTVWYBJXZ?-";
	return map[c];
}


char AASite::mapCharToNum(char c)
{
	char d;

	switch(c)
	{
		case 'A': // unambiguous
			d = 0;
			break;

		case 'C':
			d = 1;
			break;

		case 'D':
			d = 2;
			break;

		case 'E':
			d = 3;
			break;

		case 'F':
			d = 4;
			break;

		case 'G':
			d = 5;
			break;

		case 'H':
			d = 6;
			break;

		case 'I':
			d = 7;
			break;

		case 'K':
			d = 8;
			break;

		case 'L':
			d = 9;
			break;

		case 'M':
			d = 10;
			break;

		case 'N':
			d = 11;
			break;

		case 'P':
			d = 12;
			break;

		case 'Q':
			d = 13;
			break;

		case 'R':
			d = 14;
			break;

		case 'S':
			d = 15;
			break;

		case 'T':
			d = 16;
			break;

		case 'V':
			d = 17;
			break;

		case 'W':
			d = 18;
			break;

		case 'Y':
			d = 19;
			break;

		case 'B': // ambiguous
			d = 20;
			break;

		case 'J':
			d = 21;
			break;

		case 'X':
			d = 22;
			break;

		case 'Z':
			d = 23;
			break;

		case '?':
			d = 24;
			break;

		default: // missing
			d = 25;
			break;
	}

	return d;
}
