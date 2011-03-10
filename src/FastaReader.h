#ifndef FASTAREADER_H_
#define FASTAREADER_H_

#include <iostream>
#include <fstream>
#include "AlignmentReader.h"


class FastaReader: public AlignmentReader
{
public:
	FastaReader(string fileName);
	virtual ~FastaReader();
	vector<Sequence> getSequences();

private:
	ifstream _fileReader;
	string _lastLine;
};

#endif /* FASTAREADER_H_ */
