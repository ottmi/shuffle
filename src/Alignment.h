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
	Alignment();
	Alignment(int dataType);
	Alignment(Options *options);
	virtual ~Alignment();
	void addSequence(Sequence sequence);
	void removeDuplicates();
	void removeInformativeSitesDuplicates();
	void removeIncompatiblesIterative(Options *options);
	void collectSites(Options *options);
	void checkIdenticalSites();
	void computeContextIndependentScores();
	void computeContextDependentScores(unsigned int randomizations, bool writeRandomizedCo);
	void writeRandomizedCo(string prefix);
	void writeSummary(string prefix);
	Alignment getFilteredAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy);
	Alignment getInformativeSitesAlignment();
	Alignment getSubAlignment(vector<Site*> sites);
	void write(string baseName, int format);

	vector<Sequence>& getAlignment() { return _alignment; };
	Sequence getSequence(int col) { return _alignment[col]; };
	unsigned int getNumOfRows() { return _alignment.size(); };
	unsigned int getNumOfCols() { return _cols; };

private:
	void computeCo(unsigned int start, unsigned int stop, unsigned int n);
	void computePOC(unsigned int start, unsigned int stop, unsigned int n, unsigned int randomizations, bool writeRandomizedCo);
	void computeR(unsigned int start, unsigned int stop, unsigned int n);

	int _dataType;
	int _format;
	unsigned int _cols;
	vector<Sequence> _alignment;
	vector<Site*> _sites;
	vector<Site*> _informativeSites;
};

#endif /* ALIGNMENT_H_ */
