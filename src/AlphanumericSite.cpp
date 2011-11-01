#include "AlphanumericSite.h"

AlphanumericSite::AlphanumericSite(vector<Sequence>* alignment, int offset, Options *options)
{
	_unambiguousThreshold = 35;
	_type = _ALPHANUM_DATA;
	for (unsigned int i = 0; i < options->grouping.size(); i++)
		_cols.push_back(options->grouping[i] + options->groupLength * offset);
	initialize(alignment, options);
}


AlphanumericSite::AlphanumericSite(vector<unsigned int> site)
{
	_unambiguousThreshold = 35;
	_type = _ALPHANUM_DATA;
	_cols = vector<int> (1, -1);
	_site = site;
}


AlphanumericSite::~AlphanumericSite()
{
}


string AlphanumericSite::mapNumToChar(unsigned int n)
{
	string map = _ALPHANUM_MAP;
	string s;
	for (unsigned int i = 0; i < _cols.size(); i++)
	{
		s = map[n & 255] + s;
		n = n >> 8;
	}

	return s;
}


unsigned int AlphanumericSite::mapCharToNum(string s)
{
	unsigned int d = 0;
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
