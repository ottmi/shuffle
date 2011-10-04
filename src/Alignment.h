#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_
#include <vector>
#include <string>
#include "globals.h"
#include "Sequence.h"
#include "Site.h"

//using namespace std;

class Site;
class Alignment
{
public:
	Alignment(int dataType);
	Alignment(Options *options);
	virtual ~Alignment();
	void addSequence(Sequence sequence);
	void removeDuplicates();
	void collectSites(Options *options);
	void collectInformativeSites(Options *options);
	void testSymmetry(string prefix, bool extended, int windowSize, int windowStep);
	void computeCompatibilityScores(int randomizations);
	void writeSummary(string prefix);
	Alignment getModifiedAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy);
	void write(string fileName);

	vector<Sequence>* getAlignment() { return &_alignment; };
	Sequence getSequence(int col) { return _alignment[col]; };
	unsigned int getNumOfRows() { return _alignment.size(); };
	unsigned int getNumOfCols() { return _alignment[0].getLength(); };

private:
	int _dataType;
	vector<Sequence> _alignment;
	vector<Site*> _sites;
	vector<Site*> _informativeSites;
};

#endif /* ALIGNMENT_H_ */
