CXX = g++
CXXFLAGS  = -Wall -g -O3
LDFLAGS = -lm 

OBJ = AASite.o Alignment.o AlignmentReader.o DNASite.o FastaReader.o PhylipReader.o Sequence.o Site.o shuffle.o
BIN = shuffle 

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -fr $(BIN) $(OBJ)
