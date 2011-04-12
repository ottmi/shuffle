#include "globals.h"
#include "AlphanumericSite.h"

AlphanumericSite::AlphanumericSite(vector<Sequence>* alignment, vector<int> grouping, int offset)
{
	_unambiguousThreshold = 35;
	_type = _ALPHANUM_DATA;
	for (unsigned int i = 0; i < grouping.size(); i++)
		_cols.push_back(grouping[i] + offset);
	initialize(alignment);
}


AlphanumericSite::AlphanumericSite(vector<int> site)
{
	_unambiguousThreshold = 35;
	_type = _ALPHANUM_DATA;
	_cols = vector<int> (1, -1);
	_site = site;
}


AlphanumericSite::~AlphanumericSite()
{
}


string AlphanumericSite::mapNumToChar(int n)
{
	//            01234567890123456789012345678901234567
	string map = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789?-";
	string s;
	for (unsigned int i = 0; i < _cols.size(); i++)
	{
		s = map[n & 255] + s;
		n = n >> 8;
	}

	return s;
}


int AlphanumericSite::mapCharToNum(string s)
{
	int d = 0;
	for (unsigned int i = 0; i < s.size(); i++)
	{
		char c = s[i];
		d = d << 8;

		if (c >= 65 && c <= 90) // A-Z
			d+= c-65;
		else if (c >= 48 && c <= 57) // 0-9
			d+= c-48+26;
		else if (c == '?')
			d+= 36;
		else
			d+= 37;
	}
	return d;
}
