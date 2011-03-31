#include <iostream>
#include <fstream>
#include "FastaReader.h"

FastaReader::FastaReader(string fileName)
{
	cout << "FastaReader(" << fileName << ")" << endl;
	_fileReader.open(fileName.c_str());
	if (!_fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName);

	getline(_fileReader, _lastLine);
	while ((!_fileReader.eof()) && _lastLine[0] != '>')
		getline(_fileReader, _lastLine);
}


FastaReader::~FastaReader()
{
	if (!_fileReader.is_open())
		_fileReader.close();
}


vector<Sequence> FastaReader::getSequences()
{
	vector<Sequence> sequences;

	while (!_fileReader.eof())
	{
		string header;
		string seq;
		header = _lastLine;
		_lastLine = "";
		while (!_fileReader.eof() && _lastLine[0] != '>')
		{
			getline(_fileReader, _lastLine);
			if (_lastLine[0] != '>')
				seq += _lastLine;
		}
		seq = adjustString(seq);
		if (header.length() > 1 && seq.length())
		{
			Sequence s(header.substr(1), seq);
			sequences.push_back(s);
		}
	}
	return sequences;
}
