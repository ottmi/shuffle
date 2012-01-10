CXX = mpic++
CXXFLAGS  = -Wall -g -O3 -D_MPI -D_MPFR -D_DEBUG
LDFLAGS = -lm -lmpfr

OBJ = AASite.o Alignment.o AlignmentReader.o AlphanumericSite.o DNASite.o Sequence.o Site.o helper.o shuffle.o
BIN = shuffle 

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(BIN) $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(BIN) $(OBJ)
