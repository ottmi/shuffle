#ifndef PHYLIPREADER_H_
#define PHYLIPREADER_H_

#include <iostream>
#include <fstream>
#include "AlignmentReader.h"

class PhylipReader: public AlignmentReader
{
public:
	PhylipReader(string fileName);
	vector<Sequence> getSequences();

private:
	ifstream _fileReader;
	int _rows;
	int _cols;
};

#endif /* PHYLIPREADER_H_ */
