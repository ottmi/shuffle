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
};

#endif /* ALIGNMENTREADER_H_ */
