#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
using namespace std;

class Sequence {
public:
	Sequence(string name, string seq);
	virtual ~Sequence();
	string getName();
	string getSequence();
	size_t getLength();

private:
	string _name;
	string _sequence;
};

#endif /* SEQUENCE_H_ */
