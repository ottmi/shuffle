CXX = g++
CXXFLAGS  = -Wall -g -O3 -fopenmp -D_GMP
LDFLAGS = -lm -lgmpxx -lgmp -lmpfr

OBJ = AASite.o Alignment.o AlignmentReader.o AlphanumericSite.o DNASite.o FastaReader.o PhylipReader.o Sequence.o Site.o helper.o shuffle.o
BIN = shuffle 

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(BIN) $(OBJ)
