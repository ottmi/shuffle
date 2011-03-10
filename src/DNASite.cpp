#include "DNASite.h"

DNASite::DNASite(Alignment *alignment, int col)
{
/*
	 http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
	 A --> adenosine           M --> A C (amino)
	 C --> cytidine            S --> G C (strong)
	 G --> guanine             W --> A T (weak)
	 T --> thymidine           B --> G T C
	 U --> uridine             D --> G A T
	 R --> G A (purine)        H --> A C T
	 Y --> T C (pyrimidine)    V --> G C A
	 K --> G T (keto)          N --> A G C T (any)
*/
	_unambiguousCharacters = "ACGTU";
	_ambiguousCharacters = "RYKMSWBDHVN?";
	_missingCharacters = '-';

	_type = 0;
	_col = col;

	initialize(alignment->getAlignment());
}

DNASite::DNASite(string site, int col)
{
	_unambiguousCharacters = "ACGTU";
	_ambiguousCharacters = "RYKMSWBDHVN?";
	_missingCharacters = '-';
	_type = 0;
	_site = site;
	_col = col;
}


DNASite::~DNASite()
{
}
