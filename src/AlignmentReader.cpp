#include <sstream>
#include "globals.h"
#include "helper.h"
#include "AlignmentReader.h"

AlignmentReader::AlignmentReader(string fileName)
{
	_fileReader.open(fileName.c_str());
	if (!_fileReader.is_open()) throw("Error, cannot open file " + fileName);

	cout << "Reading alignment file: " << fileName << endl;

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
			cout << "The file appears to be in be in Phylip format (" << _rows << " rows, " << _cols << " columns)." << endl;
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
				stringstream s;
				s << "Unable to detect alignment format.\n" << PROGNAME << " only supports the Fasta and sequential Phylip formats.";
				throw(s.str());
			}
		}
	}
}

AlignmentReader::~AlignmentReader()
{
	if (!_fileReader.is_open()) _fileReader.close();
}

vector<Sequence> AlignmentReader::getSequences(int from, int to)
{
	string whiteSpace = " \n\t";
	vector<Sequence> sequences;

	if (from == -1) from = 1;

	if (_format == _FASTA_FORMAT)
	{
		while ((!_fileReader.eof()) && _lastLine[0] != '>')
			safeGetline(_fileReader, _lastLine);

		while (!_fileReader.eof())
		{
			string name;
			string seq;
			name = _lastLine;
			_lastLine = "";
			while (!_fileReader.eof() && _lastLine[0] != '>')
			{
				safeGetline(_fileReader, _lastLine);
				if (_lastLine[0] != '>') seq += _lastLine;
			}
			seq = adjustString(seq, false);
			name = adjustString(name.substr(1), false);
			if (name.length() > 1 && seq.length())
			{
				if (_cols && seq.length() != (unsigned int) _cols)
				{
					stringstream ss;
					ss << "Sequence #" << sequences.size() + 1 << " (" << name << ") consists of " << seq.length() << " characters. All previous sequences consisted of "
							<< _cols << " characters.";
					throw(ss.str());
				} else
					_cols = (int) seq.length();

				if (to == -1)
					sequences.push_back(Sequence(name, seq.substr(from - 1)));
				else
					sequences.push_back(Sequence(name, seq.substr(from - 1, to - from + 1)));
				_rows++;
			}
		}
	} else if (_format == _PHYLIP_FORMAT)
	{
		while (!_fileReader.eof())
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
				} else
				{
					name = _lastLine.substr(0, n);
					n = _lastLine.find_first_not_of(whiteSpace, n);
					seq = _lastLine.substr(n);
				}

				seq = adjustString(seq, false);
				if (name.length() && seq.length())
				{
					if (seq.length() != (unsigned int) _cols)
					{
						stringstream ss;
						ss << "Sequence #" << sequences.size() + 1 << " (" << name << ") consists of " << seq.length() << " characters when it should be " << _cols << ".";
						throw(ss.str());
					}

					if (to == -1)
						sequences.push_back(Sequence(name, seq.substr(from - 1)));
					else
						sequences.push_back(Sequence(name, seq.substr(from - 1, to - from + 1)));
				}
			}
		}
		if ((int) sequences.size() < _rows) cerr << "The alignment contains only " << sequences.size() << " rows, but it should be " << _rows << "." << endl;
	}

	_cols = sequences[0].getLength();

	return sequences;
}
