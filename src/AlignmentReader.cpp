#include <iostream>
#include <fstream>
#include "AlignmentReader.h"

AlignmentReader::AlignmentReader()
{
}

AlignmentReader::~AlignmentReader()
{
}

string AlignmentReader::adjustString(string s)
{
	string r = "";

	for (unsigned int i=0; i<s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != '\r' && c != ' ')
	  {
			r+= toupper(c);
	  }
	}

	return(r);
}

vector<Sequence> AlignmentReader::getSequences()
{
	vector<Sequence> dummy;
	return dummy;
}
