#include "DNASite.h"

DNASite::DNASite(vector<Sequence>* alignment, int col)
{
/*
	_unambiguousCharacters = "ACGTU";
	_ambiguousCharacters = "RYKMSWBDHVN?";
	_missingCharacters = '-';
*/

	_unambiguousThreshold = 5;
	_type = 0;
	_col = col;
	initialize(alignment);
}


DNASite::DNASite(vector<char> site, int col)
{
	_unambiguousThreshold = 5;
	_type = 0;
	_col = col;
	_site = site;
}


DNASite::~DNASite()
{
}


char DNASite::mapNumToChar(char c)
{
	string map = "ACGTURYKMSWBDHVN?-";
	return map[c];
}


char DNASite::mapCharToNum(char c)
{
	char d;

	switch(c)
	{
		case 'A':
			d = 0;
			break;

		case 'C':
			d = 1;
			break;

		case 'G':
			d = 2;
			break;

		case 'T':
			d = 3;
			break;

		case 'U':
			d = 4;
			break;

		case 'R':
			d = 5;
			break;

		case 'Y':
			d = 6;
			break;

		case 'K':
			d = 7;
			break;

		case 'M':
			d = 8;
			break;

		case 'S':
			d = 9;
			break;

		case 'W':
			d = 10;
			break;

		case 'B':
			d = 11;
			break;

		case 'D':
			d = 12;
			break;

		case 'H':
			d = 13;
			break;

		case 'V':
			d = 14;
			break;

		case 'N':
			d = 15;
			break;

		case '?':
			d = 16;
			break;

		default:
			d = 17;
			break;
	}

	return d;
}
