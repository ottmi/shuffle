#include "DNASite.h"

DNASite::DNASite(int offset, Options *options)
{
	/*
	 _unambiguousCharacters = "ACGTU";
	 _ambiguousCharacters = "RYKMSWBDHVN?";
	 _missingCharacters = '-';
	 */

	_unambiguousThreshold = 4;
	_type = 0;
	for (unsigned int i = 0; i < options->grouping.size(); i++)
		_cols.push_back(options->grouping[i] + options->groupLength * offset);
}


DNASite::DNASite(int col, vector<unsigned int> site)
{
	_unambiguousThreshold = 4;
	_type = 0;
	_cols = vector<int> (1, col);
	_site = site;
	initialize();
}


DNASite::~DNASite()
{
}


string DNASite::mapNumToChar(unsigned int n)
{
	string map = _DNA_MAP;
	string s;
	for (unsigned int i = 0; i < _cols.size(); i++)
	{
		s = map[n & 255] + s;
		n = n >> 8;
	}

	return s;
}


unsigned int DNASite::mapCharToNum(string s)
{
	unsigned int d = 0;
	for (unsigned int i = 0; i < s.size(); i++)
	{
		char c = s[i];
		d = d << 8;
		switch (c)
		{
			case 'A':
			case 'a':
				d+= 0x00;
				break;

			case 'C':
			case 'c':
				d+= 0x01;
				break;

			case 'G':
			case 'g':
				d+= 0x02;
				break;

			case 'T':
			case 't':
				d+= 0x03;
				break;

			case 'U':
			case 'u':
				d+= 0x04;
				break;

			case 'R':
			case 'r':
				d+= 0x05;
				break;

			case 'Y':
			case 'y':
				d+= 0x06;
				break;

			case 'K':
			case 'k':
				d+= 0x07;
				break;

			case 'M':
			case 'm':
				d+= 0x08;
				break;

			case 'S':
			case 's':
				d+= 0x09;
				break;

			case 'W':
			case 'w':
				d+= 0x0a;
				break;

			case 'B':
			case 'b':
				d+= 0x0b;
				break;

			case 'D':
			case 'd':
				d+= 0x0c;
				break;

			case 'H':
			case 'h':
				d+= 0x0d;
				break;

			case 'V':
			case 'v':
				d+= 0x0e;
				break;

			case 'N':
			case 'n':
				d+= 0x0f;
				break;

			case '?':
				d+= 0x10;
				break;

			default:
				d+= 0x11;
				break;
		}
	}
	return d;
}
