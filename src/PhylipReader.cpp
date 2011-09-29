#include <stdlib.h>
#include "PhylipReader.h"

string whiteSpace = " \n\t";


PhylipReader::PhylipReader(string fileName)
{
	cout << "PhylipReader(" << fileName << ")" << endl;
	_fileReader.open(fileName.c_str());
	if (! _fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName );

	string str;
	safeGetline(_fileReader, str);

	str = str.substr(str.find_first_not_of(whiteSpace));
	string rows = str.substr(0, str.find_first_of(whiteSpace));

	str = str.substr(str.find_first_of(whiteSpace));
	str = str.substr(str.find_first_not_of(whiteSpace));
	string cols = str.substr(0, str.find_first_of(whiteSpace));

	_cols = atoi(cols.c_str());
	_rows = atoi(rows.c_str());
}


vector<Sequence> PhylipReader::getSequences()
{
	vector<Sequence> sequences;
	string str;

	while (! _fileReader.eof())
    {
   		safeGetline(_fileReader, str);
   		if (str.length())
   		{
   			int n = str.find_first_of(whiteSpace);
   			string name, seq;
   			if (n == -1) // there's no whitespace, so the sequence starts at pos 11
   			{
   	   			name = str.substr(0, 10);
   	   			seq = str.substr(10);
   			}
   			else
   			{
   				name = str.substr(0, n);
   				n = str.find_first_not_of(whiteSpace, n);
   	   			seq = str.substr(n);
   			}

   			if ((int) seq.length() != _cols)
   				cerr << "Sequence #" << sequences.size() + 1 << " (" << name << ") consists of " << seq.length() << " characters when it should be " << _cols << "." << endl;
   			seq = adjustString(seq);
   			if (name.length() && seq.length())
   			{
   				Sequence s(name, seq);
   				sequences.push_back(s);
   			}
   		}
    }
	if ((int) sequences.size() < _rows)
		cerr << "The alignment has only " << sequences.size() << " rows." << endl;
    return sequences;
}
