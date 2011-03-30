BIN = shuffle 

all: $(BIN)

$(BIN):
	cd src; $(MAKE) $(MFLAGS)
	mv src/$(BIN) .

clean :
	$(RM) $(BIN)
	cd src; $(MAKE) clean
