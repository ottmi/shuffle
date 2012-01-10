#ifndef ALIGNMENTREADER_H_
#define ALIGNMENTREADER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include "Sequence.h"

using namespace std;

class AlignmentReader
{
public:
	AlignmentReader(string fileName);
	~AlignmentReader();
	vector<Sequence> getSequences();
	int getFormat() { return _format; };
	unsigned int getRows() { return _rows; };
	unsigned int getCols() { return _cols; };

private:
	ifstream _fileReader;
	string _lastLine;
	int _format; // 0=fasta, 1=phylip
	int _rows;
	int _cols;

};

#endif /* ALIGNMENTREADER_H_ */
