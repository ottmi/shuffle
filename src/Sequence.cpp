#include <stdio.h>
#include <ctype.h>
#include "Sequence.h"

string adjustString(string s)
{
	string r = "";

	for (unsigned int i=0; i<s.length(); i++)
	{
		char c = s[i];
		if (c != '\t' && c != '\n' && c != ' ')
	    {
			r+= toupper(c);
	    }
	}

	return(r);
}


Sequence::Sequence(string name, string seq)
{
	_name = name;
	_sequence = adjustString(seq);
}

Sequence::~Sequence()
{
	// TODO Auto-generated destructor stub
}

string Sequence::getName()
{
	return _name;
}

string Sequence::getSequence()
{
	return _sequence;
}

size_t Sequence::getLength()
{
	return _sequence.length();
}
