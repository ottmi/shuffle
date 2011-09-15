#include "globals.h"
#include "DNASite.h"

DNASite::DNASite(vector<Sequence>* alignment, vector<int> grouping, int offset)
{
	/*
	 _unambiguousCharacters = "ACGTU";
	 _ambiguousCharacters = "RYKMSWBDHVN?";
	 _missingCharacters = '-';
	 */

	_unambiguousThreshold = 4;
	_type = 0;
	for (unsigned int i = 0; i < grouping.size(); i++)
		_cols.push_back(grouping[i] + offset);
	initialize(alignment);
}


DNASite::DNASite(vector<int> site)
{
	_unambiguousThreshold = 5;
	_type = 0;
	_cols = vector<int> (1, -1);
	_site = site;
}


DNASite::~DNASite()
{
}


string DNASite::mapNumToChar(int n)
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


int DNASite::mapCharToNum(string s)
{
	int d = 0;
	for (unsigned int i = 0; i < s.size(); i++)
	{
		int c = s[i];
		d = d << 8;
		switch (c)
		{
			case 'A':
				d+= 0x00;
				break;

			case 'C':
				d+= 0x01;
				break;

			case 'G':
				d+= 0x02;
				break;

			case 'T':
				d+= 0x03;
				break;

			case 'U':
				d+= 0x04;
				break;

			case 'R':
				d+= 0x05;
				break;

			case 'Y':
				d+= 0x06;
				break;

			case 'K':
				d+= 0x07;
				break;

			case 'M':
				d+= 0x08;
				break;

			case 'S':
				d+= 0x09;
				break;

			case 'W':
				d+= 0x0a;
				break;

			case 'B':
				d+= 0x0b;
				break;

			case 'D':
				d+= 0x0c;
				break;

			case 'H':
				d+= 0x0d;
				break;

			case 'V':
				d+= 0x0e;
				break;

			case 'N':
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
