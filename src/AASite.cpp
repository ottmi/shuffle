#include "AASite.h"

AASite::AASite(vector<Sequence>* alignment, int offset, Options *options)
{
	/*
	 _unambiguousCharacters = "ACDEFGHIKLMNPQRSTVWY";
	 _ambiguousCharacters = "BJXZ?";
	 _missingCharacters = '-';
	 _unambiguousThreshold = 19;
	 */
	_unambiguousThreshold = 19;
	_type = 1;
	for (unsigned int i = 0; i < options->grouping.size(); i++)
		_cols.push_back(options->grouping[i] + options->groupLength * offset);

	initialize(alignment, options);
}


AASite::AASite(vector<unsigned int> site)
{
	_unambiguousThreshold = 19;
	_type = 1;
	_cols = vector<int> (1, -1);
	_site = site;
}


AASite::~AASite()
{
}


string AASite::mapNumToChar(unsigned int n)
{
	string map = _AA_MAP;
	string s;

	for (unsigned int i = 0; i < _cols.size(); i++)
	{
		s = map[n & 255] + s;
		n = n >> 8;
	}

	return s;
}


unsigned int AASite::mapCharToNum(string s)
{
	unsigned int d = 0;
	for (unsigned int i = 0; i < s.size(); i++)
	{
		char c = s[i];
		d = d << 8;
		switch (c)
		{
			case 'A': // unambiguous
				d = 0x00;
				break;

			case 'C':
				d = 0x01;
				break;

			case 'D':
				d = 0x02;
				break;

			case 'E':
				d = 0x03;
				break;

			case 'F':
				d = 0x04;
				break;

			case 'G':
				d = 0x05;
				break;

			case 'H':
				d = 0x06;
				break;

			case 'I':
				d = 0x07;
				break;

			case 'K':
				d = 0x08;
				break;

			case 'L':
				d = 0x09;
				break;

			case 'M':
				d = 0x0A;
				break;

			case 'N':
				d = 0x0B;
				break;

			case 'P':
				d = 0x0C;
				break;

			case 'Q':
				d = 0x0D;
				break;

			case 'R':
				d = 0x0E;
				break;

			case 'S':
				d = 0x0F;
				break;

			case 'T':
				d = 0x10;
				break;

			case 'V':
				d = 0x11;
				break;

			case 'W':
				d = 0x12;
				break;

			case 'Y':
				d = 0x13;
				break;

			case 'B': // ambiguous
				d = 0x14;
				break;

			case 'J':
				d = 0x15;
				break;

			case 'X':
				d = 0x16;
				break;

			case 'Z':
				d = 0x17;
				break;

			case '?':
				d = 0x18;
				break;

			default: // missing
				d = 0x19;
				break;
		}
	}

	return d;
}
