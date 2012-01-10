#include "AASite.h"

AASite::AASite(int offset, Options *options)
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
}


AASite::AASite(vector<unsigned int> site)
{
	_unambiguousThreshold = 19;
	_type = 1;
	_cols = vector<int> (1, -1);
	_site = site;
	initialize();
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
			case 'a':
				d += 0x00;
				break;

			case 'C':
			case 'c':
				d += 0x01;
				break;

			case 'D':
			case 'd':
				d += 0x02;
				break;

			case 'E':
			case 'e':
				d += 0x03;
				break;

			case 'F':
			case 'f':
				d += 0x04;
				break;

			case 'G':
			case 'g':
				d += 0x05;
				break;

			case 'H':
			case 'h':
				d += 0x06;
				break;

			case 'I':
			case 'i':
				d += 0x07;
				break;

			case 'K':
			case 'k':
				d += 0x08;
				break;

			case 'L':
			case 'l':
				d += 0x09;
				break;

			case 'M':
			case 'm':
				d += 0x0A;
				break;

			case 'N':
			case 'n':
				d += 0x0B;
				break;

			case 'P':
			case 'p':
				d += 0x0C;
				break;

			case 'Q':
			case 'q':
				d += 0x0D;
				break;

			case 'R':
			case 'r':
				d += 0x0E;
				break;

			case 'S':
			case 's':
				d += 0x0F;
				break;

			case 'T':
			case 't':
				d += 0x10;
				break;

			case 'V':
			case 'v':
				d += 0x11;
				break;

			case 'W':
			case 'w':
				d += 0x12;
				break;

			case 'Y':
			case 'y':
				d += 0x13;
				break;

			case 'B': // ambiguous
			case 'b':
				d += 0x14;
				break;

			case 'J':
			case 'j':
				d += 0x15;
				break;

			case 'X':
			case 'x':
				d += 0x16;
				break;

			case 'Z':
			case 'z':
				d += 0x17;
				break;

			case '?':
				d += 0x18;
				break;

			default: // missing
				d += 0x19;
				break;
		}
	}

	return d;
}
