#include <stdio.h>
#include <ctype.h>
#include "Sequence.h"

Sequence::Sequence(string name, string seq)
{
	_name = name;
	_sequence = seq;
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

string Sequence::getColumns(vector<int> cols)
{
	string s;

	for (unsigned int i=0; i<cols.size(); i++)
		s+= _sequence[cols[i]];

	return s;
}

size_t Sequence::getLength()
{
	return _sequence.length();
}
