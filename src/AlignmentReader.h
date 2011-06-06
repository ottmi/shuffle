#ifndef ALIGNMENTREADER_H_
#define ALIGNMENTREADER_H_

#include <vector>
#include "Sequence.h"

using namespace std;

class AlignmentReader
{
public:
	AlignmentReader();
	virtual ~AlignmentReader();
	virtual vector<Sequence> getSequences();
	string adjustString(string s, bool upercase=false);
};

#endif /* ALIGNMENTREADER_H_ */
