CXX = g++
CXXFLAGS = -Wall -g -O3
LDFLAGS = -lm

# Enable OpenMP parallelization
CXXFLAGS += -fopenmp

# Enable support for mulitple precision arithmetic
CXXFLAGS += -D_GMP
LDFLAGS += -lgmpxx -lgmp -lmpfr

OBJ = AASite.o Alignment.o AlignmentReader.o AlphanumericSite.o DNASite.o Sequence.o Site.o helper.o shuffle.o
BIN = shuffle 

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(BIN) $(OBJ)
