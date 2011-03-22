#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "FastaReader.h"
#include "PhylipReader.h"
#include "AASite.h"
#include "DNASite.h"
#include "Alignment.h"


Alignment::Alignment(int dataType)
{
	_dataType = dataType;
}


Alignment::Alignment(string fileName, int dataType)
{
	AlignmentReader *alignmentReader;

	string ext = fileName.substr(fileName.find_last_of('.') + 1);
	if (!ext.compare("phy") || !ext.compare("phylip"))
	{
		alignmentReader = new PhylipReader(fileName);
	}
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
	{
		alignmentReader = new FastaReader(fileName);
	}
	else
	{
		cerr << "Unknown input alignment format" << endl;
		exit(255);
	}
	_alignment = alignmentReader->getSequences();
	_dataType = dataType;
	delete alignmentReader;

	Site *s;
	for (unsigned int i = 0; i < getNumOfCols(); i++)
	{
		if (_dataType == 0)
			s = new DNASite(this, i);
		else
			s = new AASite(this, i);
		if (s->isInformative())
			_informativeSites.push_back(s);
		//delete s;
	}

	cout << "Alignment contains " << getNumOfRows() << " sequences with " << getNumOfCols() << " sites, " << _informativeSites.size()
			<< " of which are informative." << endl;
}

Alignment::~Alignment()
{
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
		delete _informativeSites[i];
}

void Alignment::addSequence(Sequence s)
{
	_alignment.push_back(s);
}

void Alignment::computeCompatibilityScores(int randomizations)
{
	cout << "Computing compatibility scores, doing " << randomizations << " randomizations..." << endl;

	long t1 = time(NULL);
	unsigned int total = _informativeSites.size() * randomizations;
	unsigned int count = 0;
	int myTid = 0;

#ifdef _OPENMP
	#pragma omp parallel shared(count) private(myTid)
	{
		myTid = omp_get_thread_num();
		#pragma omp for
#endif
		for (unsigned int i = 0; i < _informativeSites.size(); i++)
		{
			for (unsigned int j = i + 1; j < _informativeSites.size(); j++)
				if (_informativeSites[i]->checkCompatibility(_informativeSites[j]))
				{
					_informativeSites[i]->incComp();
					_informativeSites[j]->incComp();
				}
		}

#ifdef _OPENMP
		#pragma omp for
#endif
		for (unsigned int i = 0; i < _informativeSites.size(); i++)
		{
			_informativeSites[i]->computeCompScore(_informativeSites.size());

			int poc = 0;
			for (int r = 0; r < randomizations; r++)
			{
				int comp = 0;
				Site* randomSite = _informativeSites[i]->randomize();
				for (unsigned int j = 0; j < _informativeSites.size(); j++)
				{
					if (i != j && randomSite->checkCompatibility(_informativeSites[j]))
						comp++;
				}
				delete randomSite;
				if (_informativeSites[i]->getComp() <= comp)
					poc++;
			}

			count += randomizations;

			if (myTid == 0)
				cout << "\r" << count * 100 / total << "%" << flush;
			_informativeSites[i]->setPOC((double) poc / randomizations);
		}
#ifdef _OPENMP
	}
#endif
	long t2 = time(NULL);
	cout << "\rFinished computing scores, taking " << t2 - t1 << "s." << endl;
}

Alignment Alignment::getModifiedAlignment(double minCo, double minPOC, int maxSmin, double maxEntropy)
{
	Alignment a(_dataType);
	cout << "Creating new alignment with minCo=" << minCo << " minPOC=" << minPOC << " maxSmin=" << maxSmin << " maxEntropy=" << maxEntropy << endl;
	vector<int> sites;

	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site site = *_informativeSites[i];
		if (site.getCo() >= minCo && site.getPOC() >= minPOC && site.getSmin() <= maxSmin && site.getEntropy() <= maxEntropy)
			sites.push_back(i);
	}

	for (unsigned int i = 0; i < _alignment.size(); i++)
	{
		string seq, newSeq;
		seq = _alignment[i].getSequence();
		for (unsigned int j = 0; j < sites.size(); j++)
		{
			Site site = *_informativeSites[sites[j]];
			newSeq += seq[site.getCol()];
		}
		if (newSeq.length())
		{
			Sequence s(_alignment[i].getName(), newSeq);
			a.addSequence(s);
		}
	}
	cout << "New alignment contains " << a.getNumOfRows() << " sequences with " << a.getNumOfCols() << " sites." << endl;

	return a;
}

void Alignment::writeSummary(string fileName)
{
	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("\n\nError, cannot open file " + fileName);
	cout << "Writing site summary to " << fileName << endl;

	file << "Site No.,Smin,Entropy,OV,Co,Poc" << endl;
	for (unsigned int i = 0; i < _informativeSites.size(); i++)
	{
		Site* s = _informativeSites[i];
		file << s->getCol() + 1 << "," << s->getSmin() << "," << s->getEntropy()  << "," << s->getOV() << "," << s->getCo() << "," << s->getPOC()<< endl;
	}
}

void Alignment::write(string fileName)
{
	int type;
	string ext = fileName.substr(fileName.find_first_of('.') + 1);
	if (!ext.compare("phy") || !ext.compare("phylip"))
		type = 0;
	else if (!ext.compare("fsa") || !ext.compare("fasta"))
		type = 1;
	else
	{
		cerr << "Couldn't recognize output format." << endl;
		cerr << "Please use extension .phy or .phylip for Phylip format," << endl;
		cerr << "and .fsa or .fasta for Fasta format." << endl;
		return;
	}

	ofstream file(fileName.c_str(), ifstream::trunc);
	if (!file.is_open())
		throw("\n\nError, cannot open file " + fileName);

	switch (type)
	{
		case 0:
			cout << "Writing Phylip alignment to " << fileName << endl;
			file << getNumOfRows() << " " << getNumOfCols() << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
				file << _alignment[i].getName() << " " << _alignment[i].getSequence() << endl;
			break;

		case 1:
			cout << "Writing Fasta alignment to " << fileName << endl;
			for (unsigned int i = 0; i < getNumOfRows(); i++)
			{
				file << "<" << _alignment[i].getName() << endl;
				file << _alignment[i].getSequence() << endl;
			}
			break;
	}
	file.close();
}
