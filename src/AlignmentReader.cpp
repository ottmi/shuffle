#include <sstream>
#include <stdlib.h>
#include "globals.h"
#include "helper.h"
#include "AlignmentReader.h"

AlignmentReader::AlignmentReader(string fileName)
{
	cout << "Opening alignment file: " << fileName << endl;
	_fileReader.open(fileName.c_str());
	if (! _fileReader.is_open())
		throw("\n\nError, cannot open file " + fileName );

	_rows = 0;
	_cols = 0;

	safeGetline(_fileReader, _lastLine);

	if (_lastLine[0] == '>')
	{
		cout << "The file appears to be in Fasta format." << endl;
		_format = _FASTA_FORMAT;
	} else
	{
		stringstream ss(_lastLine);
		ss >> _rows >> _cols;
		if (_rows && _cols)
		{
			cout << "The file appears to be in be in Phylip format (" << _rows << " rows, " << _cols << " columns)."<< endl;
			_format = _PHYLIP_FORMAT;
		} else
		{
			string ext = fileName.substr(fileName.find_last_of('.') + 1);
			if (!ext.compare("fsa") || !ext.compare("fst") || !ext.compare("fasta"))
			{
				cout << "According to its extension, this file should be in Fasta format." << endl;
				_format = _FASTA_FORMAT;
			} else if (!ext.compare("phy") || !ext.compare("phylip"))
			{
				cout << "According to its extension, this file should be in Phylip format." << endl;
				_format = _PHYLIP_FORMAT;
			} else
			{
				cout << "Unable to detect alignment format." << endl;
				cout << PROGNAME << " only supports the Fasta and sequential Phylip formats."<< endl;
				exit(255);
			}
		}
	}
}


AlignmentReader::~AlignmentReader()
{
	if (!_fileReader.is_open())
		_fileReader.close();
}


vector<Sequence> AlignmentReader::getSequences()
{
	string whiteSpace = " \n\t";
	vector<Sequence> sequences;

	if (_format == _FASTA_FORMAT)
	{
		while ((!_fileReader.eof()) && _lastLine[0] != '>')
			safeGetline(_fileReader, _lastLine);

		while (!_fileReader.eof())
		{
			string header;
			string seq;
			header = _lastLine;
			_lastLine = "";
			while (!_fileReader.eof() && _lastLine[0] != '>')
			{
				safeGetline(_fileReader, _lastLine);
				if (_lastLine[0] != '>')
					seq += _lastLine;
			}
			seq = adjustString(seq, false);
			header = adjustString(header.substr(1), false);
			if (header.length() > 1 && seq.length())
			{
				Sequence s(header, seq);
				sequences.push_back(s);
				_rows++;
				if (seq.length() > (unsigned int) _cols)
					_cols = seq.length();
			}
		}
	} else if (_format == _PHYLIP_FORMAT)
	{
		while (! _fileReader.eof())
	    {
	   		safeGetline(_fileReader, _lastLine);
	   		if (_lastLine.length())
	   		{
	   			int n = _lastLine.find_first_of(whiteSpace);
	   			string name, seq;
	   			if (n == -1) // there's no whitespace, so the sequence starts at pos 11
	   			{
	   	   			name = _lastLine.substr(0, 10);
	   	   			seq = _lastLine.substr(10);
	   			}
	   			else
	   			{
	   				name = _lastLine.substr(0, n);
	   				n = _lastLine.find_first_not_of(whiteSpace, n);
	   	   			seq = _lastLine.substr(n);
	   			}

	   			if ((int) seq.length() != _cols)
	   				cerr << "Sequence #" << sequences.size() + 1 << " (" << name << ") consists of " << seq.length() << " characters when it should be " << _cols << "." << endl;
	   			seq = adjustString(seq, false);
	   			if (name.length() && seq.length())
	   			{
	   				Sequence s(name, seq);
	   				sequences.push_back(s);
	   			}
	   		}
	    }
		if ((int) sequences.size() < _rows)
			cerr << "The alignment contains only " << sequences.size() << " rows, but it should be " << _rows << "."<< endl;
	}

    return sequences;
}
